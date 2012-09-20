makeDb <- function(dbfile, textfile, tablename, sep = "\t", cutoff = 5){

	cat(file=dbfile) #create empty file
	column.names = as.character(as.matrix(read.table(textfile,sep=sep,nrows=1,header=F))) #get column names
	print(column.names)
	# create the sql statement to be used in creating database:
	for(i in 2:length(column.names)){ 
	  if(i==2) where.statement = paste(gsub("\\.","_",column.names[i]),">",cutoff,"OR")
	  else if(i==length(column.names)) where.statement = paste(where.statement,
	  	gsub("\\.","_",column.names[i]),">",cutoff)
  	  else where.statement = paste(where.statement,gsub("\\.","_",column.names[i]),">",cutoff,"OR")
  	  } #(note that "." is not an acceptable character in sql column names - is automatically replaced with "_" - so we do the same in our sql statement)
	tablename.statement = paste("main",tablename,sep=".")
	sql.statement = paste("create table",tablename.statement,"as select * from file where",where.statement)
	print(where.statement)
	print(sql.statement)
	read.csv.sql(textfile, sql=sql.statement, dbname=dbfile,sep=sep) #this is where all the action is - this creates the database, this is the function in our modified sqldf script.
	message(paste("Wrote database file",dbfile,"containing table",tablename)) # so that users know this function did something and created a file in their system.
}
