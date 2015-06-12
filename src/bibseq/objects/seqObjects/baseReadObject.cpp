#include "baseReadObject.hpp"
namespace bibseq {
baseReadObject::baseReadObject(const seqInfo& seqBase) : seqBase_(seqBase) {
  //std::cout << "baseReadObject constructor: " << std::endl;
  //std::cout << seqBase_.name_ << std::endl;
  //std::cout << seqBase_.cnt_ << std::endl;
  //std::cout << seqBase_.frac_ << std::endl;
}
void baseReadObject::printDescription(std::ostream& out, bool deep) const {
  seqBase_.printDescription(out, deep);
  out << "baseReadObject{" << std::endl << "}" << std::endl;
}
}  // namespace bib
