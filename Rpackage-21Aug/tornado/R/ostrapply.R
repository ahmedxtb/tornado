ostrapply <-
function (X, pattern, FUN = function(x, ...) x, ignore.case = FALSE, 
    ..., empty = NULL, simplify = FALSE, USE.NAMES = FALSE, combine = c) 
{
    here <- environment()
    combine <- match.funfn(combine)
    if (is.character(FUN)) {
        FUN.orig <- FUN
        FUN <- function(...) FUN.orig
    }
    else if (is.list(FUN)) {
        values.replacement <- FUN
        names.replacement <- names(FUN)
        here$FUN <- function(...) {
            idx <- match(..1, names.replacement, nomatch = match("", 
                names.replacement, nomatch = 0))
            if (idx > 0) 
                values.replacement[[idx]]
            else ..1
        }
    }
    p <- if (is.proto(FUN)) {
        FUN$X <- X
        FUN$pattern <- pattern
        FUN$simplify <- simplify
        FUN$USE.NAMES <- USE.NAMES
        FUN$combine <- combine
        proto(pre = function(this) {
            this$first <- TRUE
            this$v <- NULL
            if (!is.null(FUN[["pre"]])) 
                FUN$pre()
        }, fun = function(this, ...) {
            FUN$count <- this$count
            this$v <- if (this$first) 
                combine(FUN$fun(...))
            else c(this$v, combine(FUN$fun(...)))
            this$first <- FALSE
        }, post = function(this) {
            if (this$first) 
                this$v <- NULL
            if (!is.null(FUN[["post"]])) 
                FUN$post()
        }, )
    }
    else {
        FUN <- match.funfn(FUN)
        proto(pre = function(this) {
            this$first <- TRUE
            this$v <- NULL
        }, fun = function(this, ...) {
            this$v <- if (this$first) 
                combine(FUN(...))
            else c(this$v, combine(FUN(...)))
            this$first <- FALSE
        }, post = function(this) {
            if (this$first) 
                this$v <- NULL
        })
    }
    ff <- function(x, ...) {
        gsubfn(pattern, p, x, ignore.case = ignore.case, 
            ...)
        p$v
    }
    result <- sapply(X, ff, ..., simplify = isTRUE(simplify), 
        USE.NAMES = USE.NAMES)
    if (is.logical(simplify)) 
        result
    else {
        do.call(match.funfn(simplify), result)
    }
}
