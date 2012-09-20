gsubfn <-
function (pattern, replacement, x, backref, USE.NAMES = FALSE, 
    ignore.case = FALSE, env = parent.frame(), ...) 
{
	engine <- "R"
    R.engine <- identical(engine, "R")
    here <- environment()
    if (missing(replacement)) 
        here$replacement <- function(...) eval(parse(text = paste(..., 
            sep = "")), env)
    if (is.character(replacement)) {
        if (R.engine) 
            return(base::gsub(pattern, replacement, x, ...))
        else { stop("trying to use tcltk.  don't do that.")        }
    }
    if (is.list(replacement)) {
        values.replacement <- replacement
        names.replacement <- names(replacement)
        here$replacement <- function(...) {
            idx <- match(..1, names.replacement, nomatch = match("", 
                names.replacement, nomatch = 0))
            if (idx > 0) 
                values.replacement[[idx]]
            else ..1
        }
    }
    if (missing(pattern)) 
        pattern <- "[$]([[:alpha:]][[:alnum:].]*)|`([^`]+)`"
    pattern <- as.character(pattern)
    e <- NULL
    if (!inherits(replacement, "formula") && !is.function(replacement)) {
        e <- replacement
        e$pattern <- pattern
        e$x <- x
        e$backref <- if (!missing(backref)) 
            backref
        e$USE.NAMES <- USE.NAMES
        e$env <- env
        dots <- list(...)
        if (!is.null(names(dots))) {
            nam <- names(dots)
            for (n in nam[nam != ""]) assign(n, dots[[n]], e)
        }
        e$replacement <- function(this, ...) {
            this$count <- this$count + 1
            this$match <- c(...)
            this$fun(...)
        }
        here$replacement <- e$replacement
    }
    here$replacement <- match.funfn(replacement)
    if (missing(backref) || is.null(backref)) {
        noparen <- base::gsub("\\\\.", "", pattern)
        noparen <- base::gsub("\\[[^\\]]*\\]", "", noparen)
        backref <- -nchar(base::gsub("[^(]", "", noparen))
    }
    if (names(formals(here$replacement))[[1]] == "&") {
        backref <- abs(backref)
        if (!is.null(e)) 
            e$backref <- backref
    }
    j <- (identical(engine, "R") && !is.null(backref) && backref >= 
        0) + abs(backref)
    i <- if (!R.engine && backref >= 0) 
        0
    else 1
    j <- max(i, j)
    stopifnot(is.character(pattern), is.character(x), is.function(replacement))
    force(env)
    gsub.function <- function(x) {
        if (R.engine && !is.null(backref) && backref >= 0) {
            pattern <- paste("(", pattern, ")", sep = "")
        }
        if (!is.null(e)) {
            e$count <- 0
            if ("pre" %in% ls(e)) 
                e$pre()
        }
        repl <- function(i, j) {
            rs <- paste("\\", seq(i, j), collapse = "\002", sep = "")
            rs <- paste("\001\002", rs, "\001", sep = "")
            if (R.engine) 
                tryCatch(base::gsub(pattern, rs, x, ignore.case = ignore.case, 
                  ...), error = function(x) if (j > i) 
                  repl(i, j - 1)
                else stop(x))
            else { stop("trying to use tcltk. don't do that (2)")}
        }
        z <- repl(i, j)
        z <- strsplit(z, "\001")[[1]]
        f <- function(s) {
            if (nchar(s) > 0 && substring(s, 1, 1) == "\002") {
                s <- sub("\002$", "\002\002", s)
                L <- as.list(strsplit(s, "\002")[[1]][-1])
                do.call(replacement, L)
            }
            else s
        }
        z <- paste(sapply(z, f), collapse = "")
        if (!is.null(e) && "post" %in% ls(e)) 
            e$post()
        z
    }
    sapply(x, gsub.function, USE.NAMES = USE.NAMES)
}
