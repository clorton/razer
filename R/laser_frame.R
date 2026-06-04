# =============================================================================
# laser_frame.R — hand-written R glue layered on top of the extendr-generated
# `LaserFrame` wrapper (R/extendr-wrappers.R, which is auto-generated — do not
# edit). This file adds two ergonomics: void methods return invisibly, and
# `frame$prop` reads a property directly.
#
# Orientation for readers coming from C / C++ / C# / Python:
#   * `<-` is assignment (`=` also works but `<-` is idiomatic). Everything is a
#     value; there are no statements-vs-expressions in the C sense.
#   * Functions are first-class values: `f <- function(x) x + 1`. Assigning to a
#     specially named symbol like `` `$.LaserFrame` `` registers an S3 *method*
#     (R's ad-hoc, name-convention-based dispatch: `generic.class`).
#   * `.Call(symbol, ...)` is the C-FFI entry point — it invokes a compiled
#     routine (here, the Rust functions extendr exported) by its registered
#     symbol. The `wrap__LaserFrame__*` symbols come from extendr-wrappers.R.
#   * `x[["name"]]` extracts one element by name (no partial matching); `$` is
#     similar but does partial matching and is the operator we override below.
#   * Environments are mutable, reference-semantics namespaces (like a dict that
#     also defines variable scope). `environment(fn) <- e` rebinds the closure's
#     enclosing scope — this is how method dispatch injects `self`.
#   * `invisible(x)` returns `x` but suppresses auto-printing at the REPL (cf. a
#     `void`-ish return that still carries a value for chaining).
#   * `NULL` is the empty/absent value (closest to null / None).
# =============================================================================

# `.onLoad` is a package hook R calls automatically when the package is loaded
# (like a static initializer / module __init__). Signature is fixed by R.
#
# Wrap void-returning LaserFrame methods so they return invisibly.
#
# extendr 0.9 does not support per-method `invisible` annotations, and the
# `$.LaserFrame` dispatch mechanism replaces environment(func) on every call,
# which breaks closure-based wrappers. Instead we rewrite each function body
# in-place via body<-: the original `.Call(...)` becomes `invisible(.Call(...))`.
# The function identity (and its environment slot) is preserved.
.onLoad <- function(libname, pkgname) {
  # `c(...)` builds a vector (R's basic 1-D container); here a character vector
  # of the method names whose return value should be hidden.
  void_methods <- c(
    "add_scalar_property",
    "add_vector_property",
    "set",
    "set_col",
    "squash",
    "sort_by"
  )
  # `for (m in v)` iterates elements of the vector (Python-style foreach).
  for (m in void_methods) {
    fn <- LaserFrame[[m]]                  # pull the generated method (a function value)
    # `body(fn)` is the function's parsed body (an AST/"language" object). `call(f, args)`
    # constructs an unevaluated call node; here we wrap the old body in `invisible(...)`.
    # `body(fn) <- ...` is a replacement function — `body<-` rewrites the AST in place.
    body(fn) <- call("invisible", body(fn))
    LaserFrame[[m]] <- fn                  # store the rewritten function back on the object
  }
}

# Override `$.LaserFrame` to add direct property access: frame$prop returns
# values just like frame$get("prop"), and frame$mat_prop returns the full
# matrix just like frame$get_matrix("mat_prop").
#
# This is an S3 method: defining `` `$.LaserFrame` `` makes R dispatch here
# whenever `$` is applied to an object whose class is "LaserFrame". `self` is the
# object, `name` the (unquoted) field after the `$`, passed as a string.
#
# Method names always take priority over property names. If a property name
# collides with a method name, use $get() / $get_matrix() explicitly.
`$.LaserFrame` <- function(self, name) {
  # `return(x)` exits early with x (R also returns the last expression implicitly,
  # but explicit `return` is clearer for these guard clauses).
  # 1. Read-only frame metadata exposed as plain values, not callables
  if (name == "count")    return(.Call(wrap__LaserFrame__count,    self))
  if (name == "capacity") return(.Call(wrap__LaserFrame__capacity, self))
  # 2. Method dispatch (existing behaviour from extendr-wrappers.R). Look up a
  #    method by name; `[[` returns NULL if absent. `!is.null(fn)` == "found it".
  fn <- LaserFrame[[name]]
  if (!is.null(fn)) {
    # Rebind the method's closure environment to *this* call's environment, which
    # holds `self`. extendr's generated methods read `self` from their enclosing
    # scope, so this injection is what binds the call to this particular frame.
    environment(fn) <- environment()
    return(fn)
  }
  # 3. Scalar property shortcut. `%in%` is vectorized membership test (Python `in`).
  if (name %in% .Call(wrap__LaserFrame__scalar_names, self)) {
    return(.Call(wrap__LaserFrame__get, self, name))
  }
  # 4. Vector property shortcut — returns (count × ncols) matrix
  if (name %in% .Call(wrap__LaserFrame__vector_names, self)) {
    return(.Call(wrap__LaserFrame__get_matrix, self, name))
  }
  NULL    # not metadata, method, or property -> absent (last expr is the return value)
}

# Make `frame[["x"]]` behave exactly like `frame$x` by pointing the `[[` S3
# method at the same function value (functions are first-class, so this is a
# plain alias — no wrapper). The roxygen tags above each export route the docs
# to the shared LaserFrame help topic; `@usage NULL` suppresses an auto-usage line.
#' @rdname LaserFrame
#' @usage NULL
#' @export
`[[.LaserFrame` <- `$.LaserFrame`

# Assignment form: `frame$prop <- value` dispatches here (note the `$<-` name).
# An S3 replacement method must take `value` as its last argument and return the
# (possibly modified) object — R rebinds the variable to whatever this returns.
#' @rdname LaserFrame
#' @usage NULL
#' @export
`$<-.LaserFrame` <- function(self, name, value) {
  if (name %in% .Call(wrap__LaserFrame__scalar_names, self)) {
    .Call(wrap__LaserFrame__set, self, name, value)   # write through to the Rust column
    # Return self invisibly. The frame is an external pointer with reference
    # semantics, so the write already mutated it; we return self to satisfy the
    # `$<-` contract without printing it.
    return(invisible(self))
  }
  if (name %in% .Call(wrap__LaserFrame__vector_names, self)) {
    # `stop(...)` raises an R error (like `throw`). `paste0` concatenates strings
    # with no separator (`paste` uses a space); `\"` escapes a literal quote.
    stop(paste0(
      "direct assignment to vector property '", name,
      "' is not supported; use $set_col(\"", name, "\", col, values)"
    ))
  }
  stop(paste0("'", name, "' is not a property of this LaserFrame"))
}
