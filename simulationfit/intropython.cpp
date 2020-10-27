#include "func.hpp"
#include <boost/python.hpp>

using namespace std;

class MyClass {
public:
        void add(int valueA, int valueB) {
                _result = valueA + valueB;
        }
        int value() {
                return _result;
        }
private:
        int _result;
};

BOOST_PYTHON_MODULE(MyClassModule) {
        using namespace boost::python;

        class_<MyClass>("MyClass")
                .def("add", &MyClass::add)
                .def("value", &MyClass::value);
}
