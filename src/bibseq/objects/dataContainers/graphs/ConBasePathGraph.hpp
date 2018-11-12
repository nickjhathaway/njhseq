#pragma once
/*
 * ConBasePathGraph.hpp
 *
 *  Created on: Jan 14, 2017
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

#include "njhseq/utils.h"

namespace njhseq {






class ConBasePathGraph {
public:
	class ConPath {
	public:
		struct PosBase {
			uint32_t pos_;
			char base_;
			std::string getUid() const;
			Json::Value toJson() const;
		};

		ConPath(double count);

		double count_ = 1;
		std::vector<PosBase> bases_;
		void addPosBase(uint32_t pos, char base);
		void addPosBase(const PosBase & pb);

		std::string getUid() const;

		Json::Value toJson() const;
	};

	class edge;
	class node {
	public:

		node(const ConPath::PosBase & val,
				double cnt, double frac) ;
		ConPath::PosBase val_;
		double cnt_;
		double frac_;

		std::vector<std::shared_ptr<edge>> headEdges_;
		std::vector<std::shared_ptr<edge>> tailEdges_;

		uint32_t visitCount_ = 0;

		void resetVisitCount();

		void addHead(const std::shared_ptr<edge> & e);

		void addTail(const std::shared_ptr<edge> & e) ;

		bool headless() const;

		bool tailless() const;

		void addToWritingPath(std::ostream & out,
				std::string currentPath);
		void addToPath(std::vector<ConPath> & paths,
				ConPath currentPath);
		Json::Value toJson() const;
	};
	class edge {
	public:
		edge(const std::shared_ptr<node> & head,
				const std::shared_ptr<node> & tail,
				double cnt);
		std::weak_ptr<node> head_;
		std::weak_ptr<node> tail_;
		double cnt_;
		Json::Value toJson() const;
	};

	std::unordered_map<std::string, std::shared_ptr<node>> nodes_;

	void addNode(const ConPath::PosBase & n, double cnt, double frac) ;

	void addEdge(const std::string & head,
			const std::string & tailName,
			double cnt) ;

	void writePaths(std::ostream & out) const;

	std::vector<ConPath> getPaths() const ;

	Json::Value createSankeyOutput() const ;
};




}  // namespace njhseq
