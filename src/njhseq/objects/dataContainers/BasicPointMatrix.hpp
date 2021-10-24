#pragma once
/*
 * BasicPointMatrix.hpp
 *
 *  Created on: Jun 13, 2017
 *      Author: nick
 */
// njhseq - A library for analyzing sequence data
// Copyright (C) 2012-2018 Nicholas Hathaway <nicholas.hathaway@umassmed.edu>,
//
// This file is part of njhseq.
//
// njhseq is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// njhseq is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with njhseq.  If not, see <http://www.gnu.org/licenses/>.
//
#include "njhseq/objects/dataContainers/graphs/UndirWeightedGraph.hpp"
#include "njhseq/concurrency/PairwisePairFactory.hpp"

namespace njhseq {


template<typename VAL = double>
class BasicPointMatrix {
public:
	/**@brief a class wrapping a vector to represent a basic point
	 *
	 */
	class BasicPoint {
	public:
		/**@brief the values of the points
		 *
		 * @param values the values to store in the point
		 */
		BasicPoint(const std::vector<VAL> & values) :
				vals_(values) {
		}
		std::vector<VAL> vals_; /**< the values for this point */
		/**@brief get the euclidean distance between two points
		 *
		 * @param p1 point 1
		 * @param p2 point 2
		 * @return the distance between the two points
		 */
		static double euDist(const BasicPoint& p1, const BasicPoint & p2) {
			//dangerous, no guarantee that these will be the same length;
			double sum = 0;
			for (const auto pos : iter::range(p1.vals_.size())) {
				sum += std::pow(p1.vals_[pos] - p2.vals_[pos], 2);
			}
			return std::sqrt(sum);
		}
		/**@brief get the distacne between this point and another point
		 *
		 * @param p2 the other point
		 * @return the distance
		 */
		double euDist(const BasicPoint & p2) const {
			return euDist(*this, p2);
		}
		/**@brief write out the values stored in the point (doesn't add newlines or and adds to current line)
		 *
		 * @param out the out stream to write
		 * @param delim the delimiter to put inbetween the points
		 */
		void writeVals(std::ostream & out, const std::string & delim = "\t") const {
			out << njh::conToStr(vals_, delim);
		}
	};
	/**@brief a simple matrix class, mostly invented to do dbscan on a matrix of points
	 *
	 * @param dbPars the db scans parameters to run latter
	 */
	BasicPointMatrix(
			const typename njhUndirWeightedGraph<double,
					std::shared_ptr<BasicPointMatrix<VAL>::BasicPoint>>::dbscanPars & dbPars) :
			dbscanPars_(dbPars) {
	}
	const typename njhUndirWeightedGraph<double,
			std::shared_ptr<BasicPointMatrix<VAL>::BasicPoint>>::dbscanPars dbscanPars_; /**< the dbscan parameters to run dbscan latter */
	std::vector<std::shared_ptr<BasicPoint>> points_; /**< the matrix */
	VecStr rowNames_; /**< possible row names, could be empty if no row names */
	VecStr colNames_; /**< possible column names, could be empty if no column names */
	std::unique_ptr<
			njhUndirWeightedGraph<double,
					std::shared_ptr<BasicPointMatrix<VAL>::BasicPoint>>> graph_; /**< the unweighed graph to create from the matrix */
	/**@brief set the unweighed graph from the input matrix and dbscan parameters
	 *
	 * @param numThreads number of threads to use
	 * @param verbose whether to be verbose when running
	 */
	void setGraph(uint32_t numThreads, bool verbose) {
		if (0 == numThreads) {
			numThreads = 1;
		}

		graph_ = std::make_unique<
				njhUndirWeightedGraph<double,
						std::shared_ptr<BasicPointMatrix<VAL>::BasicPoint>>>();

		for (const auto & pos : iter::range(points_.size())) {
			graph_->addNode(estd::to_string(pos), points_[pos]);
		}

		uint32_t belowEp = 0;
		/**@todo this appears to be actually fairly slow, i think it's mostly because the eu calculations is so fast, perhaps a better way of multithreading this can be done
		 *
		 */
		PairwisePairFactory pairFactory(points_.size());
		uint32_t pairBatchCount = 100000;
		std::mutex graphMut;
		struct PairDist {
			PairDist(const PairwisePairFactory::PairwisePair & pair, double dist) :
					pair_(pair), dist_(dist) {
			}
			PairwisePairFactory::PairwisePair pair_;
			double dist_;
		};

		std::function<void()> addToGraph =
				[&graphMut, &pairFactory,&pairBatchCount,&belowEp,this]() {
					PairwisePairFactory::PairwisePairVec pairs;
					std::vector<PairDist> belowEps;
					while(pairFactory.setNextPairs(pairs, pairBatchCount)) {
						for(const auto & pair : pairs.pairs_) {
							auto dist = points_[pair.row_]->euDist(*points_[pair.col_]);
							if (dist < dbscanPars_.eps_) {
								belowEps.emplace_back(PairDist{pair, dist});
							}
						}
					}
					if(!belowEps.empty()) {
						std::lock_guard<std::mutex> lock(graphMut);
						belowEp += belowEps.size();
						for(const auto & bEps : belowEps) {
							graph_->addEdge(estd::to_string(bEps.pair_.row_), estd::to_string(bEps.pair_.col_),
									bEps.dist_);
						}
					}
				};

		njh::concurrent::runVoidFunctionThreaded(addToGraph, numThreads);

		if (verbose) {
			std::cout << std::endl;
			std::cout << "below: " << belowEp << "/" << pairFactory.totalCompares_ << std::endl;
		}
	}
	/**@brief a check to throw if the graph hasn't been set yet
	 *
	 */
	void checkGraphThrow() const {
		if (nullptr == graph_) {
			std::stringstream ss;
			ss << __PRETTY_FUNCTION__
					<< ", error trying to write graph when it isn't set" << "\n";
			throw std::runtime_error { ss.str() };
		}
	}
	/**@brief write out the graph to a file defined by opts
	 *
	 * @param opts the out options
	 */
	void writeGraph(const OutOptions & opts) const {
		checkGraphThrow();
		std::ofstream outFile;
		opts.openFile(outFile);
		writeGraph(outFile);
	}
	/**@brief write graph out to a stream
	 *
	 * write noise points as -2
	 *
	 * @param out the out stream to write to
	 */
	void writeGraph(std::ostream & out) const {
		checkGraphThrow();
		for (const auto & nPos : iter::range(graph_->nodes_.size())) {
			if (!rowNames_.empty()) {
				out << rowNames_[nPos] << "\t";
			}
			const auto & n = graph_->nodes_[nPos];
			n->value_->writeVals(out, "\t");
			if (std::numeric_limits<uint32_t>::max() == n->group_) {
				//noise group
				out << "\t" << -2;
			} else {
				out << "\t" << n->group_;
			}
			out << std::endl;
		}
	}
	/**@brief a factory to create a matrix to do dbscan analysis on from a file
	 *
	 * @param fnp the matrix file to read from, can have row names (first row is 1 less than the rest of the rows)
	 * @param dbPars the parameters for db scan
	 * @return a basic matrix
	 */
	static BasicPointMatrix<VAL> readInBasicMatrix(const bfs::path & fnp,
			const typename njhUndirWeightedGraph<double,
					std::shared_ptr<BasicPointMatrix<VAL>::BasicPoint>>::dbscanPars & dbPars,
					bool firstColRowNames = false,
					bool colNames = false) {
		BasicPointMatrix<VAL> ret(dbPars);
		if (!bfs::exists(fnp)) {
			std::stringstream ss;
			ss << __PRETTY_FUNCTION__ << ":" << fnp << " doesn't exist" << "\n";
			throw std::runtime_error { ss.str() };
		}
		std::ifstream inFile(fnp.string());
		if (!inFile) {
			std::stringstream ss;
			ss << __PRETTY_FUNCTION__ << ": error in opening " << fnp << "\n";
			throw std::runtime_error { ss.str() };
		}
		bool inputMightHaveRowNames = njh::files::hasPossibleRowNames(fnp);
		uint32_t expectedColNum = njh::files::getExpectedNumCol(fnp);
		std::string line;
		if(inputMightHaveRowNames){
			colNames = true;
			firstColRowNames = true;
		}
		if(colNames){
			njh::files::crossPlatGetline(inFile, line);
			ret.colNames_= tokenizeString(line, "whitespace");
		}
		while (njh::files::crossPlatGetline(inFile, line)) {
			auto toks = tokenizeString(line, "whitespace");
			uint32_t rowNameOffset = (firstColRowNames ? 1 : 0);
			if (toks.size() != expectedColNum + rowNameOffset) {
				std::stringstream ss;
				ss << __PRETTY_FUNCTION__ << ": error in processing line " << line
						<< "\n";
				ss << "Was expecting " << expectedColNum + rowNameOffset
						<< " columns not " << toks.size() << "\n";
				throw std::runtime_error { ss.str() };
			}
			if (firstColRowNames) {
				ret.rowNames_.emplace_back(toks[0]);
			}
			ret.points_.emplace_back(
					std::make_shared<BasicPoint>(
							vecStrToVecNum<double>(
									VecStr { toks.begin() + rowNameOffset, toks.end() })));
		}
		return ret;
	}
};

}  // namespace njhseq


//example program below
//int runDbscan(const njh::progutils::CmdArgs & inputCommands){
//	bfs::path fnp = "";
//	OutOptions outOpts(bfs::path("out.tab.txt"));
//	njhUndirWeightedGraph<double, std::shared_ptr<BasicPointMatrix<double>::BasicPoint>>::dbscanPars dbscanPars;
//	dbscanPars.eps_ = 1;
//	dbscanPars.minEpNeighbors_ = 5;
//	uint32_t numThreads = 1;
//	seqSetUp setUp(inputCommands);
//	setUp.processVerbose();
//	setUp.setOption(fnp, "--fnp", "Name of the input matrix", true);
//	setUp.processWritingOptions(outOpts);
//	setUp.setOption(dbscanPars.eps_, "--eps", "Epsilon (distance sensitivity of algorithm)");
//	setUp.setOption(dbscanPars.minEpNeighbors_, "--minpts", "The minimum number of epsilon neighbors");
//	setUp.setOption(numThreads, "--numThreads", "Number of threads");
//	setUp.finishSetUp(std::cout);
//	std::ofstream outFile;
//	outOpts.openFile(outFile);
//	njh::stopWatch watch;
//	watch.setLapName("Reading in");
//	auto mat = BasicPointMatrix<double>::readInBasicMatrix(fnp, dbscanPars);
//	watch.startNewLap("Adding nodes");
//
//	mat.setGraph(numThreads, setUp.pars_.verbose_);
//
//	watch.startNewLap("dbscan");
//	mat.graph_->dbscan(dbscanPars);
//
//	watch.startNewLap("output");
//	mat.writeGraph(outFile);
//
//	if(setUp.pars_.verbose_){
//		watch.logLapTimes(std::cout, true, 6, true);
//	}
//	return 0;
//}



