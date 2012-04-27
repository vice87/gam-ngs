#ifndef _MERGE_IN_CUT_TAIL_FAILED_HPP_
#define _MERGE_IN_CUT_TAIL_FAILED_HPP_

#include <string>
#include <exception>

class MergeInCutTailFailed : public std::exception
{
    std::string _what;

  public:
    MergeInCutTailFailed(const std::string& what):
                                 _what(what) {}

    MergeInCutTailFailed(const char* what):
                                 _what(std::string(what)) {}

    ~MergeInCutTailFailed() throw() {}

    virtual const char* what() const throw() { return _what.c_str(); }
};


#endif // _MERGE_IN_CUT_TAIL_FAILED_