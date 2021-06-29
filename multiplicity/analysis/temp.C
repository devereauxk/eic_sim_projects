void temp() {
	std::string a  = "a";
	const char * b = "b";
	a.append("b\n");
	printf("%s", a.c_str());
	return 0;
}
