.onLoad <- function(libname, pkgname) {
    ExperimentHub::createHubAccessors(pkgname, "PscanR_backgrounds_v2")
}
