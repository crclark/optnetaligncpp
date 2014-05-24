#include<iostream>
#include<boost/filesystem/operations.hpp>

namespace bfs=boost::filesystem;
int main()
{
bfs::path p("second.cpp");
if(bfs::exists(p))
std::cout<<p.leaf()<<std::endl;
}