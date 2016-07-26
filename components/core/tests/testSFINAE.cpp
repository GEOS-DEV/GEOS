
#ifdef __clang__
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wglobal-constructors"
#pragma clang diagnostic ignored "-Wexit-time-destructors"
#endif

#include "gtest/gtest.h"

#ifdef __clang__
#pragma clang diagnostic push
#endif

#include "codingUtilities/sfinae.hpp"


struct Foo_MemberData
{
  int memberName = 1;
};
struct Bar_MemberData : Foo_MemberData {};

struct Foo_StaticMemberData
{
  static int memberName;
};
int Foo_StaticMemberData::memberName = 1;
struct Bar_StaticMemberData : Foo_StaticMemberData {};


struct Foo_MemberFunction_1Arg
{
  int memberName(int a) { return a; }
};
struct Bar_MemberFunction_1Arg : Foo_MemberFunction_1Arg {};

struct Foo_MemberFunction_2Arg
{
  double memberName2(int a,double b) { return a+b; }
};
struct Bar_MemberFunction_2Arg : Foo_MemberFunction_2Arg {};

struct Foo_StaticMemberFunction_1Arg
{
  static int memberName( int a) { return a; }
};
struct Bar_StaticMemberFunction_1Arg : Foo_StaticMemberFunction_1Arg {};


struct Foo_Using
{
  using memberName = int;
};
struct Bar_Using : Foo_Using {};

struct Foo_Typedef
{
  typedef int memberName;
};
struct Bar_Typedef : Foo_Typedef {};

struct Foo_EnumClass
{
  enum class memberName
  {
    enum1,
    enum2
  };
};
struct Bar_EnumClass : Foo_EnumClass {};


HAS_MEMBER_DATA(memberName)
HAS_STATIC_MEMBER_DATA(memberName)
HAS_MEMBER_FUNCTION( memberName, int, , VA_LIST(int), VA_LIST(int(1)) )
HAS_MEMBER_FUNCTION( memberName2, double, , VA_LIST(int,double), VA_LIST(int(1),double(1)) )
HAS_STATIC_MEMBER_FUNCTION(memberName,int,int(1))
HAS_ENUM(memberName)
HAS_ALIAS(memberName)

TEST(test_sfinae,test_has_datamember)
{
  EXPECT_TRUE( has_datamember_memberName<Foo_MemberData>::value );
  EXPECT_FALSE( has_datamember_memberName<Foo_StaticMemberData>::value );
  EXPECT_FALSE( has_datamember_memberName<Foo_MemberFunction_1Arg>::value );
  EXPECT_FALSE( has_datamember_memberName<Foo_StaticMemberFunction_1Arg>::value );
  EXPECT_FALSE( has_datamember_memberName<Foo_Using>::value );
  EXPECT_FALSE( has_datamember_memberName<Foo_Typedef>::value );
  EXPECT_FALSE( has_datamember_memberName<Foo_EnumClass>::value );
}

TEST(test_sfinae,test_has_staticdatamember)
{
  EXPECT_FALSE( has_staticdatamember_memberName<Foo_MemberData>::value );
  EXPECT_TRUE(  has_staticdatamember_memberName<Foo_StaticMemberData>::value );
  EXPECT_FALSE( has_staticdatamember_memberName<Foo_MemberFunction_1Arg>::value );
  EXPECT_FALSE( has_staticdatamember_memberName<Foo_StaticMemberFunction_1Arg>::value );
  EXPECT_FALSE( has_staticdatamember_memberName<Foo_Using>::value );
  EXPECT_FALSE( has_staticdatamember_memberName<Foo_Typedef>::value );
  EXPECT_FALSE( has_staticdatamember_memberName<Foo_EnumClass>::value );
}

TEST(test_sfinae,test_has_memberfunction)
{
  EXPECT_FALSE( has_memberfunction_memberName<Foo_MemberData>::value );
  EXPECT_FALSE( has_memberfunction_memberName<Foo_StaticMemberData>::value );
  EXPECT_TRUE(  has_memberfunction_memberName<Foo_MemberFunction_1Arg>::value );
  EXPECT_FALSE( has_memberfunction_memberName<Foo_StaticMemberFunction_1Arg>::value );
  EXPECT_FALSE( has_memberfunction_memberName<Foo_Using>::value );
  EXPECT_FALSE( has_memberfunction_memberName<Foo_Typedef>::value );
  EXPECT_FALSE( has_memberfunction_memberName<Foo_EnumClass>::value );

  EXPECT_TRUE(  has_memberfunction_memberName2<Foo_MemberFunction_2Arg>::value );

}

TEST(test_sfinae,test_has_staticmemberfunction)
{
  EXPECT_FALSE( has_staticmemberfunction_memberName<Foo_MemberData>::value );
  EXPECT_FALSE( has_staticmemberfunction_memberName<Foo_StaticMemberData>::value );
  EXPECT_FALSE(  has_staticmemberfunction_memberName<Foo_MemberFunction_1Arg>::value );
  EXPECT_TRUE( has_staticmemberfunction_memberName<Foo_StaticMemberFunction_1Arg>::value );
  EXPECT_FALSE( has_staticmemberfunction_memberName<Foo_Using>::value );
  EXPECT_FALSE( has_staticmemberfunction_memberName<Foo_Typedef>::value );
  EXPECT_FALSE( has_staticmemberfunction_memberName<Foo_EnumClass>::value );
}

TEST(test_sfinae,test_has_enum)
{
  EXPECT_FALSE( has_enum_memberName<Foo_MemberData>::value );
  EXPECT_FALSE( has_enum_memberName<Foo_StaticMemberData>::value );
  EXPECT_FALSE( has_enum_memberName<Foo_MemberFunction_1Arg>::value );
  EXPECT_FALSE( has_enum_memberName<Foo_StaticMemberFunction_1Arg>::value );
  EXPECT_FALSE( has_enum_memberName<Foo_Using>::value );
  EXPECT_FALSE( has_enum_memberName<Foo_Typedef>::value );
  EXPECT_TRUE(  has_enum_memberName<Foo_EnumClass>::value );
}


TEST(test_sfinae,test_has_alias)
{
  EXPECT_FALSE( has_alias_memberName<Foo_MemberData>::value );
  EXPECT_FALSE( has_alias_memberName<Foo_StaticMemberData>::value );
  EXPECT_FALSE(  has_alias_memberName<Foo_MemberFunction_1Arg>::value );
  EXPECT_FALSE( has_alias_memberName<Foo_StaticMemberFunction_1Arg>::value );
  EXPECT_TRUE(  has_alias_memberName<Foo_Using>::value );
  EXPECT_TRUE(  has_alias_memberName<Foo_Typedef>::value );
  EXPECT_FALSE( has_alias_memberName<Foo_EnumClass>::value );
}
