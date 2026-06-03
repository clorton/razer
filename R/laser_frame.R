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

# Override `$.LaserFrame` to add direct property access: frame$prop returns
# values just like frame$get("prop"), and frame$mat_prop returns the full
# matrix just like frame$get_matrix("mat_prop").
#
# Method names always take priority over property names. If a property name
# collides with a method name, use $get() / $get_matrix() explicitly.
`$.LaserFrame` <- function(self, name) {
  # 1. Read-only frame metadata exposed as plain values, not callables
  if (name == "count")    return(.Call(wrap__LaserFrame__count,    self))
  if (name == "capacity") return(.Call(wrap__LaserFrame__capacity, self))
  # 2. Method dispatch (existing behaviour from extendr-wrappers.R)
  fn <- LaserFrame[[name]]
  if (!is.null(fn)) {
    environment(fn) <- environment()
    return(fn)
  }
  # 3. Scalar property shortcut
  if (name %in% .Call(wrap__LaserFrame__scalar_names, self)) {
    return(.Call(wrap__LaserFrame__get, self, name))
  }
  # 4. Vector property shortcut — returns (count × ncols) matrix
  if (name %in% .Call(wrap__LaserFrame__vector_names, self)) {
    return(.Call(wrap__LaserFrame__get_matrix, self, name))
  }
  NULL
}

#' @rdname LaserFrame
#' @usage NULL
#' @export
`[[.LaserFrame` <- `$.LaserFrame`

#' @rdname LaserFrame
#' @usage NULL
#' @export
`$<-.LaserFrame` <- function(self, name, value) {
  if (name %in% .Call(wrap__LaserFrame__scalar_names, self)) {
    .Call(wrap__LaserFrame__set, self, name, value)
    return(invisible(self))
  }
  if (name %in% .Call(wrap__LaserFrame__vector_names, self)) {
    stop(paste0(
      "direct assignment to vector property '", name,
      "' is not supported; use $set_col(\"", name, "\", col, values)"
    ))
  }
  stop(paste0("'", name, "' is not a property of this LaserFrame"))
}
