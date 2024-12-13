library(dplyr)
library(shiny)
library(DECIPHER)
# DECIPHER::BrowseSeqs , but returns the HTML
GetBrowseHTML <- function (
    myXStringSet, 
    htmlFile = paste(tempdir(), "/myXStringSet.html", sep = ""), 
    openURL = interactive(), 
    colorPatterns = TRUE, 
    highlight = NA, 
    patterns = c("-", alphabet(myXStringSet, baseOnly = TRUE)), 
    colors = substring(rainbow(length(patterns), v = 0.8, start = 0.9, end = 0.7), 1, 7), 
    colWidth = Inf, 
    ...) 
{
  if (!is(myXStringSet, "XStringSet")) 
    stop("myXStringSet must be an XStringSet.")
  if (length(myXStringSet) == 0) 
    stop("No sequence information to display.")
  if (is(patterns, "DNAStringSet")) {
    if (!is(myXStringSet, "DNAStringSet")) 
      stop("patterns must be a DNAStringSet.")
    type <- 1L
  }
  else if (is(patterns, "RNAStringSet")) {
    if (!is(myXStringSet, "RNAStringSet")) 
      stop("patterns must be a RNAStringSet.")
    type <- 2L
  }
  else if (is(patterns, "AAStringSet")) {
    if (!is(myXStringSet, "AAStringSet")) 
      stop("patterns must be a AAStringSet.")
    type <- 3L
  }
  else if (is(patterns, "list")) {
    type <- -1L
  }
  else {
    type <- 0L
  }
  if (type > 0L) {
    patterns <- as.character(patterns)
    if (type == 1L) {
      patterns <- gsub("M", "[ACM]", patterns, fixed = TRUE)
      patterns <- gsub("R", "[AGR]", patterns, fixed = TRUE)
      patterns <- gsub("W", "[ATW]", patterns, fixed = TRUE)
      patterns <- gsub("S", "[CGS]", patterns, fixed = TRUE)
      patterns <- gsub("Y", "[CTY]", patterns, fixed = TRUE)
      patterns <- gsub("K", "[GTK]", patterns, fixed = TRUE)
      patterns <- gsub("V", "[ACGMRSV]", patterns, fixed = TRUE)
      patterns <- gsub("H", "[ACTMWYH]", patterns, fixed = TRUE)
      patterns <- gsub("D", "[AGTRWKD]", patterns, fixed = TRUE)
      patterns <- gsub("B", "[CGTSYKB]", patterns, fixed = TRUE)
      patterns <- gsub("N", "[ACGTMRWSYKVHDBN]", patterns, 
                       fixed = TRUE)
    }
    else if (type == 2L) {
      patterns <- gsub("M", "[ACM]", patterns, fixed = TRUE)
      patterns <- gsub("R", "[AGR]", patterns, fixed = TRUE)
      patterns <- gsub("W", "[AUW]", patterns, fixed = TRUE)
      patterns <- gsub("S", "[CGS]", patterns, fixed = TRUE)
      patterns <- gsub("Y", "[CUY]", patterns, fixed = TRUE)
      patterns <- gsub("K", "[GUK]", patterns, fixed = TRUE)
      patterns <- gsub("V", "[ACGMRSV]", patterns, fixed = TRUE)
      patterns <- gsub("H", "[ACUMWYH]", patterns, fixed = TRUE)
      patterns <- gsub("D", "[AGURWKD]", patterns, fixed = TRUE)
      patterns <- gsub("B", "[CGUSYKB]", patterns, fixed = TRUE)
      patterns <- gsub("N", "[ACGUMRWSYKVHDBN]", patterns, 
                       fixed = TRUE)
    }
    else {
      patterns <- gsub("B", "[NDB]", patterns, fixed = TRUE)
      patterns <- gsub("Z", "[QEZ]", patterns, fixed = TRUE)
      patterns <- gsub("J", "[ILJ]", patterns, fixed = TRUE)
      patterns <- gsub("X", "[ARNDCQEGHILKMFPSTWYVUOBJZX]", 
                       patterns, fixed = TRUE)
      patterns <- gsub("*", "\\*", patterns, fixed = TRUE)
    }
  }
  else if (type == 0) {
    if (!is.null(patterns) && !is.character(patterns)) 
      stop("patterns must be a character vector.")
    if (any(grepl("=|\"|<|>|[1-9]|[a-z]", patterns))) 
      stop("patterns cannot contain numbers, lower case characters, or the characters (=, <, >, \").")
  }
  if (any(patterns == "")) 
    stop("patterns cannot be empty (i.e., '').")
  if (type < 0) {
    if (length(myXStringSet) != length(patterns)) 
      stop("patterns is not the same length as myXStringSet.")
  }
  else {
    w <- which(patterns %in% c("?", "*", "+", "."))
    if (length(w) > 0) 
      patterns[w] <- paste("\\", patterns[w], sep = "")
    if (length(colors) != length(patterns) || !is.character(colors)) 
      stop("colors must be a character vector of the same length as patterns.")
  }
  if (is.character(htmlFile)) {
    htmlfile <- file(htmlFile, "w")
    on.exit(close(htmlfile))
  }
  else if (!inherits(htmlFile, "connection")) {
    stop("htmlFile must be a character string or connection.")
  }
  if (!is.logical(openURL) || is.na(openURL)) 
    stop("openURL must be TRUE or FALSE.")
  if (!is.logical(colorPatterns) && !is.numeric(colorPatterns)) 
    stop("colorPatterns must be a logical or numeric.")
  if (is.numeric(colorPatterns)) {
    if ((length(colorPatterns)%%2) == 1 || length(colorPatterns) == 
        0) 
      stop("colorPatterns must specify all start and endpoints.")
    if (any((colorPatterns[seq(2, length(colorPatterns), 
                               2)] - colorPatterns[seq(1, length(colorPatterns), 
                                                       2)]) < 0)) 
      stop("colorPatterns specifies a negative range.")
    if (any(colorPatterns <= 0)) 
      stop("colorPatterns must be a positive numeric.")
    if (length(colorPatterns) > 2) 
      if (any((colorPatterns[seq(3, length(colorPatterns), 
                                 2)] - colorPatterns[seq(2, length(colorPatterns) - 
                                                         2, 2)]) <= 0)) 
        stop("Ranges specified in colorPatterns must be non-overlapping.")
    if (max(colorPatterns) > max(width(myXStringSet))) 
      stop("Ranges specified in colorPatterns are out-of-bounds.")
  }
  if (is.numeric(colorPatterns) & !is.infinite(colWidth)) 
    stop("colWidth must be Inf if colorPatterns is numeric.")
  if (is.null(names(myXStringSet))) {
    names(myXStringSet) <- 1:length(myXStringSet)
  }
  else {
    names(myXStringSet) <- iconv(names(myXStringSet), from = "", 
                                 to = "UTF-8", sub = " ")
    names(myXStringSet) <- gsub("\t", " ", names(myXStringSet), 
                                fixed = TRUE)
    names(myXStringSet) <- gsub("\n", " ", names(myXStringSet), 
                                fixed = TRUE)
  }
  if (!is.na(highlight)) {
    if (highlight < 0 || highlight > length(myXStringSet)) 
      stop("highlight must be 0 or the index of a sequence in myXStringSet.")
    if (highlight != floor(highlight)) 
      stop("highlight be be a whole number.")
  }
  if (is(myXStringSet, "DNAStringSet") || is(myXStringSet, 
                                             "RNAStringSet") || is(myXStringSet, "AAStringSet")) {
    myXStringSet <- c(myXStringSet, ConsensusSequence(myXStringSet, 
                                                      ...))
  }
  else {
    myXStringSet <- c(myXStringSet, as(consensusString(myXStringSet, 
                                                       ...), class(myXStringSet)))
  }
  l <- length(myXStringSet)
  names(myXStringSet)[l] <- "Consensus"
  html <- as.character(myXStringSet)
  if (type < 0) {
    s <- strsplit(html, "", fixed = TRUE)
    n <- nchar(html)
    v <- vector("list", length(html) - 1)
    for (i in seq_len(length(html) - 1)) {
      if (!is.matrix(patterns[[i]])) 
        stop("All elements of the list patterns contain a single matrix.")
      if (nrow(patterns[[i]]) != 3) 
        stop("All elements of the list patterns must be a matrix with 3 rows.")
      if (any(patterns[[i]] < 0) || any(patterns[[i]] > 
                                        1)) 
        stop("All elements of patterns[[", i, "]] are not between 0 and 1.")
      v[[i]] <- character(n[i])
      if (any(patterns[[i]] > 1) || any(patterns[[i]] < 
                                        0)) 
        stop("All values of patterns[[", i, "]] must be between 0 and 1.")
      w <- which(s[[i]] %in% c(LETTERS, letters, "*"))
      if (length(w) != ncol(patterns[[i]])) 
        stop("The number of columns in patterns[[", i, 
             "]] is different than the number of letters in myXStringSet[", 
             i, "]")
      if (length(w) > 0) 
        v[[i]][w] <- apply(patterns[[i]], 2, function(x) {
          rgb(x[1], x[2], x[3])
        })
      if (is.numeric(colorPatterns)) {
        start <- 0L
        cPs <- c(colorPatterns, length(v[[i]]) + 1)
        for (j in seq_len(length(cPs))) {
          if (j%%2) {
            end <- cPs[j]
            if (end > length(v[[i]])) 
              end <- length(v[[i]]) + 1L
            if (start < end - 1) 
              v[[i]][(start + 1L):(end - 1L)] <- ""
          }
          else {
            start <- cPs[j]
            if (start > length(v[[i]]) + 1) 
              break
          }
        }
      }
    }
  }
  if (!is.na(highlight)) {
    if (highlight == 0) {
      highlight <- length(html)
      index <- 1:(length(html) - 1L)
    }
    else {
      index <- (1:(length(html) - 1L))[-highlight]
    }
    html <- sapply(html, strsplit, split = "", fixed = TRUE)
    for (i in index) {
      L <- min(length(html[[highlight]]), length(html[[i]]))
      w <- which(html[[i]][1:L] == html[[highlight]][1:L])
      if (length(w) > 0) 
        html[[i]][w] <- "·"
    }
    html <- sapply(html, paste, collapse = "")
  }
  maxW <- max(width(myXStringSet))
  if (maxW == 0) 
    stop("No sequence information to display.")
  if (colWidth > maxW) 
    colWidth <- maxW
  for (i in seq_len(l)) {
    html[i] <- paste(html[i], paste(rep(" ", maxW - nchar(html[i])), 
                                    collapse = ""), sep = "")
  }
  if (maxW < 20) {
    if (maxW < 10) {
      counter <- maxW
    }
    else {
      counter <- 10
    }
  }
  else {
    counter <- 20
  }
  offset <- (counter - 1) - nchar(maxW)
  if (offset < 0) 
    offset <- 0
  legend <- paste(paste(rep(" ", offset), collapse = ""), format(seq(counter, 
                                                                     maxW, by = counter)), collapse = "")
  counter <- ceiling(counter/2)
  tickmarks <- paste(paste(rep("'", counter - 1), collapse = ""), 
                     rep("|", floor(maxW/counter)), collapse = "", sep = "")
  tickmarks <- paste(tickmarks, paste(rep("'", maxW - counter * 
                                            floor(maxW/counter)), collapse = ""), sep = "")
  starts <- seq(1, maxW, colWidth)
  stops <- starts + colWidth - 1L
  stops[length(stops)] <- maxW
  temp <- character((length(html) + 5L) * length(starts))
  count <- 1L
  for (i in 1:length(starts)) {
    temp[count:(count + 1L)] <- substring(c(legend, tickmarks), 
                                          starts[i], stops[i])
    count <- count + 2L
    temp[c(count:(count + length(html) - 2L), count + length(html))] <- substring(html, 
                                                                                  starts[i], stops[i])
    count <- count + length(html) + 3L
  }
  html <- temp
  if (length(starts) > 1) {
    lengths <- numeric(length(starts) * l)
    myLengths <- character(length(starts) * (l + 5L))
    for (i in 1:length(starts)) {
      s <- substring(myXStringSet, starts[i], stops[i])
      lengths[((i - 1L) * l + 1L):(i * l)] <- nchar(s) - 
        rowSums(letterFrequency(BStringSet(s), c("-", 
                                                 " ", ".", "+")))
      if (i > 1) 
        lengths[((i - 1L) * l + 1L):(i * l)] <- lengths[((i - 
                                                            1L) * l + 1L):(i * l)] + lengths[((i - 2L) * 
                                                                                                l + 1L):((i - 1L) * l)]
      myLengths[((i - 1L) * (l + 5L) + 3L):(i * (l + 5L) - 
                                              4L)] <- lengths[((i - 1L) * l + 1L):(i * l - 
                                                                                     1L)]
      myLengths[(i * (l + 5L) - 2L)] <- lengths[i * l]
    }
  }
  else {
    lengths <- width(myXStringSet) - rowSums(letterFrequency(myXStringSet, 
                                                             c("-", ".", "+")))
    myLengths <- c("", "", lengths[-seq(l, length(lengths), 
                                        l)], "", lengths[seq(l, length(lengths), l)], "", 
                   "")
  }
  if (is.numeric(colorPatterns) || colorPatterns) {
    if (type < 0) {
      s <- strsplit(html, "", fixed = TRUE)
      for (i in seq_along(v)) {
        count <- 2L
        for (j in seq_along(starts)) {
          w <- which(v[[i]][starts[j]:stops[j]] != "")
          if (length(w) > 0) 
            s[[i + count]][w] <- paste("<span style=\"color:#FFF;background:", 
                                       v[[i]][w + starts[j] - 1], "\">", s[[i + 
                                                                              count]][w], "</span>", sep = "")
          count <- count + length(myXStringSet) + 5L
        }
      }
      html <- sapply(s, paste, collapse = "")
    }
    else {
      patterns <- paste("(", patterns, ifelse(nchar(patterns) == 
                                                1, "+)", ")"), sep = "")
      classes <- paste("<span class=\"_", seq_along(patterns), 
                       "\">\\1</span>", sep = "")
      if (is.numeric(colorPatterns)) {
        htm <- substring(html, 0, colorPatterns[1] - 
                           1)
        for (i in seq(1, length(colorPatterns), 2)) {
          htmi <- substring(html, colorPatterns[i], colorPatterns[i + 
                                                                    1])
          for (j in seq_along(colors)) {
            htmi <- gsub(patterns[j], classes[j], htmi)
          }
          end <- ifelse(i == (length(colorPatterns) - 
                                1), max(nchar(html)), colorPatterns[i + 2] - 
                          1)
          htm <- paste(htm, htmi, substring(html, colorPatterns[i + 
                                                                  1] + 1, end), sep = "")
        }
        html <- htm
      }
      else {
        for (j in seq_along(colors)) {
          html <- gsub(patterns[j], classes[j], html)
        }
      }
    }
    html <- paste(html, myLengths, "", sep = "    ")
    html <- paste("", format(c("", "", names(myXStringSet)[1:(l - 
                                                                1)], "", "Consensus", "", ""), justify = "right"), 
                  html, sep = "    ")
    html <- gsub("&", "&amp;", html, fixed = TRUE)
    styles <- character()
    for (i in seq_along(colors)) {
      styles <- paste(styles, "span._", i, " {background:", 
                      colors[i], "; color:#FFF;} ", sep = "")
    }
    styles <- paste("<style type=text/css> ", styles, "</style>", 
                    sep = "")
    html <- c("<html><head><meta http-equiv=\"Content-Type\" content=\"text/html; charset=utf-8\"></head>", 
              styles, "<pre>", html, "</pre></html>")
  }
  else {
    html <- paste(html, myLengths, "", sep = "    ")
    html <- paste("", format(c("", "", names(myXStringSet)[1:(l - 
                                                                1)], "", "Consensus", "", ""), justify = "right"), 
                  html, sep = "    ")
    html <- c("<html>", "<pre>", html, "</pre></html>")
  }
  html <- gsub("·", "&#183;", html, fixed = TRUE)
 return(html)
}
# Global variables can go here
if (!file.exists("data.RData")){
  #J01859.1
  coli <-gsub("\n", "", "AAATTGAAGAGTTTGATCATGGCTCAGATTGAACGCTGGCGGCAGGCCTAACACATGCAAGTCGAACGGT
AACAGGAAGAAGCTTGCTCTTTGCTGACGAGTGGCGGACGGGTGAGTAATGTCTGGGAAACTGCCTGATG
GAGGGGGATAACTACTGGAAACGGTAGCTAATACCGCATAACGTCGCAAGACCAAAGAGGGGGACCTTCG
GGCCTCTTGCCATCGGATGTGCCCAGATGGGATTAGCTAGTAGGTGGGGTAACGGCTCACCTAGGCGACG
ATCCCTAGCTGGTCTGAGAGGATGACCAGCCACACTGGAACTGAGACACGGTCCAGACTCCTACGGGAGG
CAGCAGTGGGGAATATTGCACAATGGGCGCAAGCCTGATGCAGCCATGCCGCGTGTATGAAGAAGGCCTT
CGGGTTGTAAAGTACTTTCAGCGGGGAGGAAGGGAGTAAAGTTAATACCTTTGCTCATTGACGTTACCCG
CAGAAGAAGCACCGGCTAACTCCGTGCCAGCAGCCGCGGTAATACGGAGGGTGCAAGCGTTAATCGGAAT
TACTGGGCGTAAAGCGCACGCAGGCGGTTTGTTAAGTCAGATGTGAAATCCCCGGGCTCAACCTGGGAAC
TGCATCTGATACTGGCAAGCTTGAGTCTCGTAGAGGGGGGTAGAATTCCAGGTGTAGCGGTGAAATGCGT
AGAGATCTGGAGGAATACCGGTGGCGAAGGCGGCCCCCTGGACGAAGACTGACGCTCAGGTGCGAAAGCG
TGGGGAGCAAACAGGATTAGATACCCTGGTAGTCCACGCCGTAAACGATGTCGACTTGGAGGTTGTGCCC
TTGAGGCGTGGCTTCCGGAGCTAACGCGTTAAGTCGACCGCCTGGGGAGTACGGCCGCAAGGTTAAAACT
CAAATGAATTGACGGGGGCCCGCACAAGCGGTGGAGCATGTGGTTTAATTCGATGCAACGCGAAGAACCT
TACCTGGTCTTGACATCCACGGAAGTTTTCAGAGATGAGAATGTGCCTTCGGGAACCGTGAGACAGGTGC
TGCATGGCTGTCGTCAGCTCGTGTTGTGAAATGTTGGGTTAAGTCCCGCAACGAGCGCAACCCTTATCCT
TTGTTGCCAGCGGTCCGGCCGGGAACTCAAAGGAGACTGCCAGTGATAAACTGGAGGAAGGTGGGGATGA
CGTCAAGTCATCATGGCCCTTACGACCAGGGCTACACACGTGCTACAATGGCGCATACAAAGAGAAGCGA
CCTCGCGAGAGCAAGCGGACCTCATAAAGTGCGTCGTAGTCCGGATTGGAGTCTGCAACTCGACTCCATG
AAGTCGGAATCGCTAGTAATCGTGGATCAGAATGCCACGGTGAATACGTTCCCGGGCCTTGTACACACCG
CCCGTCACACCATGGGAGTGGGTTGCAAAAGAAGTAGGTAGCTTAACCTTCGGGAGGGCGCTTACCACTT
TGTGATTCATGACTGGGGTGAAGTCGTAACAAGGTAACCGTAGGGGAACCTGCGGTTGGATCACCTCCTT
A")
  # region definitions: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4802574/ Table 1
  regions<- read.csv(text="Region	ignore1	ignore2	start	end
V1	8	789	8	96
V2	790	2697	97	306
V3	2698	4069	307	487
V4	4070	7044	488	746
V5	7045	9533	747	885
V6	9534	10454	886	1029
V7	10455	12258	1030	1180
V8	12259	13597	1181	1372
V9	13598	14371	1373	1468", sep="\t") |> dplyr::rowwise() |> dplyr::mutate(seq=substr(start = start, stop = end, coli))
  primers_raw <- read.csv("Abellan-Schneyder-2021-primers.tsv", sep="\t", col.names=c("Region", "Forward", "Reverse","Fseq", "Rseq", "Specificity", "Temp", "Ref"), header = FALSE)  |> 
    dplyr::mutate(Fseq=gsub(" *", "", Fseq),
                  Rseq=gsub(" *", "", Rseq))
                  
  primers_r <- primers_raw |> dplyr::select(Region, Reverse, Rseq) |> dplyr::rename(Forward=Reverse) 
  primers_r$Fseq <- as.character(Biostrings::reverseComplement(DNAStringSet(primers_r$Rseq)))
  primers <- dplyr::bind_rows(
    primers_raw, 
    primers_r
    ) |> 
    dplyr::rename("Name"="Forward", "Seq"="Fseq") |> 
    dplyr::select(Region, Name, Seq) |> 
    dplyr::mutate(Seq=gsub(" *", "", Seq),
                  Name=paste(ifelse(grepl("F$", Name), gsub("(.+?)-(.*)", "\\1", Region), gsub("(.+?)-(.+)", "\\2", Region)), Name)) |> 
    dplyr::select(Name, Seq) |> 
    dplyr::distinct()
  # get rc of rev primers
  primerseqs<- c( coli,regions$seq, primers$Seq)
  names(primerseqs) <- c( "E. coli Reference", regions$Region, primers$Name)
  primerseqs <- Biostrings::DNAStringSet(primerseqs)
  save(list = c("coli", "primerseqs"), file = "data.RData")
}else{
  load("data.RData")
}

n <- 200

#alignment <- DECIPHER::AlignSeqs(primerseqs, )
# Define the UI
ui <- bootstrapPage(
  checkboxGroupInput("refs", "References:",inline = TRUE,
                     choiceNames  = names(primerseqs),
                     choiceValues = names(primerseqs), 
                     selected = names(primerseqs)[c(1, 19, 20)],
  ),
  textInput('in1', 'Primer1', value = "AGAGTTTGATYMTGGCTCAG"),
  textInput('in2', 'Primer2', value = "AAGTCGTAACAAGGTARCCGTA"),
  textInput('rc', 'reverse complementer'),
  shiny::verbatimTextOutput('rc'),
  # verbatimTextOutput('reg'),
  shiny::uiOutput('aln')
)




# Define the server code
server <- function(input, output) {
  alignment <- reactive({
    newseqs <- c(input$in1,input$in2)
    print(input$refs)
    allseqs <- xscat(c(primerseqs[input$refs], newseqs))
    names(allseqs) <- c(names(primerseqs[input$refs]), "Primer1", "Primer2")    
    pairs = data.frame(Pattern=1, Subject=seq_along(allseqs))
    tmpaln <- DECIPHER::AlignPairs(
      gapOpening = -5,gapExtension = -10,
      allseqs, allseqs, pairs=pairs, perfectMatch = 1, misMatch = -3, type="both"
    )
    # names(tmpaln[[3]])
    # names(tmpaln) <- c(names(primerseqs[input$refs]), 
    #                     paste("Primer1", ScoreAlignment(tmpaln[c(1,which(names(tmpaln) == "Primer1"))],  substitutionMatrix=as.matrix(nucs) *1)), 
    #                     paste("Primer2", ScoreAlignment(tmpaln[c(1,which(names(tmpaln) == "Primer2"))],  substitutionMatrix=as.matrix(nucs) *1)))
    tmpaln[[3]]
  })
  thisalphabet<- alphabet(primerseqs, baseOnly = TRUE)
output$aln <- renderUI(HTML(paste(
  GetBrowseHTML(
    alignment(),
    highlight = 1,
    patterns = c("-", thisalphabet),
    colors=c("#EDEDED", "#CC7A00", "#29CC00", "#00CCCC", "#2900CC"),
    ), collapse = "\n")))
rc <- reactive({
  req(input$rc)
  as.character(Biostrings::reverseComplement(DNAStringSet(input$rc)))
})
output$rc <- renderText(rc())
output$reg <-   renderText("    E. coli Reference    ")
}

# Return a Shiny app object
shinyApp(ui = ui, server = server)

