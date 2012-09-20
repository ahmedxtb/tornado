strapply <-
function (X, pattern, FUN = function(x, ...) x, backref = NULL, 
    ..., empty = NULL, ignore.case = FALSE, perl = TRUE, engine = "R", 
    simplify = FALSE, USE.NAMES = FALSE, combine = c) 
{
    combine <- match.funfn(combine)
    stopifnot(!missing(pattern))
    pattern <- as.character(pattern)
    if (is.proto(FUN) || perl) 
        engine <- "R"
    if (identical(engine, "R")) 
        return(ostrapply(X = X, ignore.case = ignore.case, pattern = pattern, 
            FUN = FUN, backref = backref, ..., empty = empty, 
            perl = perl, simplify = simplify, USE.NAMES = USE.NAMES, 
            combine = combine))
    if (is.proto(FUN)) {
    }
    else if (is.character(FUN)) {
        FUN.orig <- FUN
        FUN <- function(...) FUN.orig
    }
    else if (is.list(FUN)) {
        values.replacement <- FUN
        names.replacement <- names(FUN)
        FUN <- function(...) {
            idx <- match(..1, names.replacement, nomatch = match("", 
                names.replacement, nomatch = 0))
            if (idx > 0) 
                values.replacement[[idx]]
            else ..1
        }
    }
    else {
        FUN <- match.funfn(FUN)
    }
    ff <- function(x) {
        s <- strapply1(x, pattern, backref, ignore.case)
        if (length(s) == 0 && !is.null(empty)) 
            s <- matrix(empty, 1)
        L <- lapply(seq_len(ncol(s)), function(j) {
            combine(do.call(FUN, as.list(s[, j])))
        })
        do.call("c", L)
    }
    result <- sapply(X, ff, simplify = is.logical(simplify) && 
        simplify, USE.NAMES = USE.NAMES)
    if (is.logical(simplify)) 
        result
    else do.call(match.funfn(simplify), result)
}
