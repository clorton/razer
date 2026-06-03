# Wrap void-returning LaserFrame methods so they return invisibly.
#
# extendr 0.9 does not support per-method `invisible` annotations, and the
# `$.LaserFrame` dispatch mechanism replaces environment(func) on every call,
# which breaks closure-based wrappers. Instead we rewrite each function body
# in-place via body<-: the original `.Call(...)` becomes `invisible(.Call(...))`.
# The function identity (and its environment slot) is preserved.
.onLoad <- function(libname, pkgname) {
  void_methods <- c(
    "add_scalar_property",
    "add_vector_property",
    "set",
    "set_col",
    "squash",
    "sort_by"
  )
  for (m in void_methods) {
    fn <- LaserFrame[[m]]
    body(fn) <- call("invisible", body(fn))
    LaserFrame[[m]] <- fn
  }
}
