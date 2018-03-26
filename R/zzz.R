.First.lib <-function(lib,pkg)
{
ver <- read.dcf(file.path(lib, pkg, "DESCRIPTION"), "Version")
     ver <- as.character(ver)	

library.dynam("nlt",pkg,lib)

cat("nlt", ver, "loaded\n")
}
