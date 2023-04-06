R__LOAD_LIBRARY(libeicsmear)

void make_tree(std::string filstr){

  erhic::DisKinematics::BoundaryWarning=false;
  
  std::string dirstr = "outfiles";
  std::string inputstr = dirstr + "/" + filstr;
  BuildTree(inputstr,dirstr);
}
