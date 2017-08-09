#include <gtest/gtest.h>
#include "dataRepository/ManagedGroup.hpp"
#include "dataRepository/ViewWrapper.hpp"
#include "common/DataTypes.hpp"

namespace geosx {
namespace dataRepository {


TEST(testSidre, simpleRestore) {
  ManagedGroup * root = new ManagedGroup(std::string("int64_data"), nullptr);
    // ViewWrapper<int64_array> data_view = root.registerViewWrapper("data");
    // data_view.resize(10);
    // int64_ptr data = data_view.data();
    
    // for (int i = 0; i < 10; i++) {
    //     data[i] = i;
    // }

  EXPECT_TRUE(false);
}

} // end namespace dataRepository
} // end namespace goesx