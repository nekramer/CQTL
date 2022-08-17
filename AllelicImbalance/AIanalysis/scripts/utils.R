# Function to parse variant IDs in the format chrX:pos:ref:alt and
# convert to GRanges

varID_to_GRanges <- function(varIDs){
  
  chr <- varIDs %>% str_split(":") %>% map(1) %>% unlist()
  pos <- varIDs %>% str_split(":") %>% map(2) %>% unlist()
  
  ranges <- GRanges(seqnames = chr, 
                    ranges = IRanges(start = as.numeric(pos), 
                                     end = as.numeric(pos)))
 return(ranges) 
}
