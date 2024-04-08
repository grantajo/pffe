#ifndef UTIL_EXCEPTIONS_H
#define UTIL_EXCEPTIONS_H
#include <exception>

namespace Util
{
namespace Exception
{
class UnitTest : public std::runtime_error
{
public:
    using std::runtime_error::runtime_error;
};

class IO : public std::runtime_error
{
public:
    using std::runtime_error::runtime_error;
};
class Runtime : public std::runtime_error
{
public:
    using std::runtime_error::runtime_error;
};


class Numeric : public std::runtime_error
{
public:
    using std::runtime_error::runtime_error;
};


}
}

#endif
