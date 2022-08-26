#include "unittests_common.hh"

#if __cplusplus >= 201103L || _MSC_VER >= 1800

#include <OpenVolumeMesh/System/MemoryInclude.hh>

TEST(MakeUniqueTest, MakeUniqueTest) {
  std::unique_ptr<int> foo;
  auto bar = ptr::make_unique<int>(5);
  foo = std::move(bar);

  EXPECT_EQ(*foo, 5);
  EXPECT_EQ(bar.get(), nullptr);
}

#endif
