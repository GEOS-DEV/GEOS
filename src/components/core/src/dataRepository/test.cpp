#include "ManagedGroup.hpp"
#inclue "ViewWrapper.hpp"

namespace geosx {
namespace dataRepository {

bool testSimpleRestore() {
    ManagedGroup root = new ManagedGroup("int64_data", nullptr);
    ViewWrapper<int64_array> data_view = root.registerViewWrapper("data");
    data_view.resize(10);
    int64_ptr data = data_view.data();
    
    for (int i = 0; i < 10; i++) {
        data[i] = i;
    }

    return true; 
}

int main(int argc, char const *argv[]) {
    bool result = testSimpleRestore();
    cout << result << endl;
    return 0;
}

} //end namespace dataRepository
} //end namespace geosx