read.csv.sql <-
function (file, sql = "select * from file", header = TRUE, sep = ",", 
    row.names, eol, skip, filter, nrows, field.types, comment.char, 
    dbname = tempfile(), drv = "SQLite", ...) 
{
	require(proto)
    file.format <- list(header = header, sep = sep)
    if (!missing(eol)) 
        file.format <- append(file.format, list(eol = eol))
    if (!missing(row.names)) 
        file.format <- append(file.format, list(row.names = row.names))
    if (!missing(skip)) 
        file.format <- append(file.format, list(skip = skip))
    if (!missing(filter)) 
        file.format <- append(file.format, list(filter = filter))
    if (!missing(nrows)) 
        file.format <- append(file.format, list(nrows = nrows))
    if (!missing(field.types)) 
        file.format <- append(file.format, list(field.types = field.types))
    if (!missing(comment.char)) 
        file.format <- append(file.format, list(comment.char = comment.char))
    pf <- parent.frame()
    if (missing(file) || is.null(file) || is.na(file)) 
        file <- ""
    tf <- NULL
    if (substring(file, 1, 7) == "http://" || substring(file, 
        1, 6) == "ftp://") {
        tf <- tempfile()
        on.exit(unlink(tf), add = TRUE)
        download.file(file, tf, mode = "wb")
        file <- tf
    }
    p <- proto(pf, file = file(file))
    p <- do.call(proto, list(pf, file = file(file)))
    sqldf(sql, envir = p, file.format = file.format, dbname = dbname, 
        drv = drv, ...)
}
