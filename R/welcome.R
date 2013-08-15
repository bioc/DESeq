# ==========================================================================
# package initialization
# ==========================================================================
.onAttach = function(libname, pkgname) {
    msg = "Welcome to 'DESeq'. For improved performance, usability and functionality, please consider migrating to 'DESeq2'."
    msg = strwrap(msg, exdent=4, indent=4)
    packageStartupMessage(paste(msg, collapse="\n"), appendLF=TRUE)
}

