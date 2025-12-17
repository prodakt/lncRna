# ===================================================================
# LncRna Package: Example Analysis Pipeline
# ===================================================================
#
# Ten skrypt demonstruje kompletny przepływ pracy analizy lncRNA przy
# użyciu funkcji z pakietu 'lncRna'.
#
# --- PREREQUISITES ---
# 1. Zainstaluj pakiet 'lncRna' i wszystkie jego zależności:
#    if (!require("BiocManager")) install.packages("BiocManager")
#    BiocManager::install("prodakt/lncRna@test")
#
# 2. Pobierz i rozpakuj dane do analizy do podkatalogu o nazwie 'data'.
#
# 3. Upewnij się, że Twoim katalogiem roboczym jest folder, w którym
#    znajduje się ten skrypt oraz podkatalog 'data'.
#
# ===================================================================


# ===================================================================
# Sekcja 1: Wczytywanie Danych Wejściowych
# ===================================================================
# W tej sekcji zdefiniujemy ścieżki do plików i wczytamy wszystkie
# niezbędne dane do środowiska R.

message("Sekcja 1: Wczytywanie danych...")
getwd()
setwd("F:/Dragoni/Prodakt/lncRna/")
# --- 1.1: Zdefiniuj ścieżki do plików ---
# Używamy `file.path` dla kompatybilności z różnymi systemami operacyjnymi.
dataDir <- "data"
stringtieGtfFile <- file.path(dataDir, "stringtie_merged.gtf.gz")
refGtfFile <- file.path(dataDir, "Mus_musculus.GRCm39.107.gtf.gz")
geneCountsFile <- file.path(dataDir, "gene_count_matrix.csv")
transcriptCountsFile <- file.path(dataDir, "transcript_count_matrix.csv")
degsFile <- file.path(dataDir, "DEGs.csv")
delsFile <- file.path(dataDir, "DELs.csv")
cdsFastaFile <- file.path(dataDir, "Mus_musculus.GRCm39.cds.all.fa.gz")
ncRnaFastaFile <- file.path(dataDir, "Mus_musculus.GRCm39.ncrna.fa.gz")

# --- 1.2: Wczytaj pliki GTF ---
stringtieGtf <- rtracklayer::import.gff(stringtieGtfFile)
refGtf <- rtracklayer::import.gff(refGtfFile)

# --- 1.3: Wczytaj macierze ekspresji ---
genesCounts <- read.csv(geneCountsFile, row.names = 1)
rownames(genesCounts) <- gsub("\\|.*", "", rownames(genesCounts)) # Oczyść nazwy genów
transcriptsCounts <-  read.table("data/transcript_count_matrix.csv", header = TRUE, sep = ",", row.names = 1)

# --- 1.4: Wczytaj listy genów DE ---
degs <- read.csv2(degsFile, row.names = 1)
dels <- read.csv2(delsFile, row.names = 1)

# --- 1.5: Wczytaj sekwencje referencyjne ---
cdsReferenceSeqs <- seqinr::read.fasta(cdsFastaFile, seqtype = "DNA", as.string = TRUE)
ncReferenceSeqs <- seqinr::read.fasta(ncRnaFastaFile, seqtype = "DNA", as.string = TRUE)


# ===================================================================
# Sekcja 2: Wstępne Filtrowanie i Przygotowanie Danych
# ===================================================================
message("Sekcja 2: Filtrowanie i przygotowanie danych...")

# --- 2.1: Wyekstrahuj statystyki transkryptów z pliku GTF StringTie ---
transcriptStats <- getGtfStats(gtfObject = stringtieGtf)

# --- 2.2: Zidentyfikuj znane biotypy z referencyjnego GTF ---
knownBiotypes <- getRefBiotypes(gtfObject = refGtf, level = "transcript")
length(knownBiotypes$transcript_biotype)
table(knownBiotypes$transcript_biotype)
knownLncRnaIds <- knownBiotypes[knownBiotypes$transcript_biotype == "lncRNA", "transcript_id"]
knownPcRnaIds <- knownBiotypes[knownBiotypes$transcript_biotype %in% "protein_coding",]$transcript_id
head(knownPcRnaIds)
length(knownPcRnaIds)
# --- 2.3: Filtruj transkrypty na podstawie cech i ekspresji ---
# Usuń znane transkrypty kodujące białka
potentialLncRnaTranscripts <- transcriptStats[!(transcriptStats$transcript_id %in% knownPcRnaIds), ]
head(potentialLncRnaTranscripts)
length(potentialLncRnaTranscripts$transcript_id)
# Filtruj po długości (>200 nt) i liczbie egzonów (>1)
potentialLncRnaTranscriptsIds <- potentialLncRnaTranscripts[potentialLncRnaTranscripts$exons > 1 & potentialLncRnaTranscripts$transLength >= 200,]$transcript_id
head(potentialLncRnaTranscriptsIds)
length(potentialLncRnaTranscriptsIds)
# Zidentyfikuj transkrypty z wystarczającą ekspresją
expressedTranscriptIds <- rownames(transcriptsCounts[rowSums(transcriptsCounts) > 10, ])
# Finalna lista potencjalnych lncRNA (po wszystkich filtrach)
potentialLncRnaIds <- potentialLncRnaTranscripts$transcript_id[potentialLncRnaTranscripts$transcript_id %in% expressedTranscriptIds]
head(potentialLncRnaIds)
length(potentialLncRnaIds)
message(paste("Znaleziono", length(potentialLncRnaIds), "potencjalnych lncRNA po filtrowaniu."))


# ===================================================================
# Sekcja 3: Przygotowanie Zestawów Treningowych i Testowych
# ===================================================================
message("Sekcja 3: Przygotowanie zestawów treningowych i testowych...")

# --- 3.1: Podziel sekwencje referencyjne na zbiory testowe i treningowe ---
set.seed(12345) # Dla odtwarzalności
cdsSets <- splitTestTrain(ids = cdsReferenceSeqs, trainProportion = 0.6)
ncSets <- splitTestTrain(ids = ncReferenceSeqs, trainProportion = 0.6)

# --- 3.2: Przygotuj sekwencje do predykcji potencjału kodującego ---
# Będziemy przewidywać dla naszych kandydatów ORAZ dla znanych sekwencji testowych
transcriptome <- seqinr::read.fasta(file.path(dataDir, "merged_transtriptome.fa.gz"), seqtype = "DNA", as.string = TRUE)
seqsToPredict <- c(
  transcriptome[names(transcriptome) %in% potentialLncRnaIds],
  cdsReferenceSeqs[names(cdsReferenceSeqs) %in% cdsSets$test],
  ncReferenceSeqs[names(ncReferenceSeqs) %in% ncSets$test]
)


# ===================================================================
# Sekcja 4: Analiza Potencjału Kodującego
# ===================================================================
message("Sekcja 4: Analiza potencjału kodującego...")

# --- 4.1: Wczytaj i połącz wyniki z narzędzi ---
codpotResults <- codPotToTbl(
  cpc2Outfile      = file.path(dataDir, "lncCodPot_MMus/CPC2_Mm_lnc.txt.txt"),
  feelncOutfile    = file.path(dataDir, "lncCodPot_MMus/FEELnc_codpot_RF.txt"),
  cnciOutfile      = file.path(dataDir, "lncCodPot_MMus/CNCI.index"),
  cpatOutfile      = file.path(dataDir, "lncCodPot_MMus/CPAT_Mm_lnc.txt"),
  cpatCutoff       = 0.64,
  plekOutfile      = file.path(dataDir, "lncCodPot_MMus/PLEK_Mm_lnc.txt"),
  lncFinderOutfile = file.path(dataDir, "lncCodPot_MMus/LncFinder_results_more5.csv")
)

# --- 4.2: Wizualizuj zgodność narzędzi (opcjonalnie) ---
if (requireNamespace("venn", quietly = TRUE)) {
  # Wykres dla wszystkich narzędzi
  plotCodPotVenn(codPot = codpotResults)
  # Wykres dla wybranych narzędzi
  plotCodPotVenn(codPot = codpotResults, selection = c(1, 0, 0, 1, 1))
}

# --- 4.3: Przygotuj dane do oceny wydajności narzędzi ---
singleToolPerformance <- sumSingleTools(
  codPotList = codpotResults,
  ncTest = ncSets$test,
  cdsTest = cdsSets$test
)

# --- 4.4: Ocena wydajności pojedynczych narzędzi ---
selectedTools <- c("CPC2", "PLEK", "FEELnc", "CPAT", "CNCI", "LncFinder")
bestToolResults <- bestTool(
  sumSingleToolsList = singleToolPerformance,
  tools = selectedTools
)
print("Wydajność pojedynczych narzędzi:")
print(bestToolResults)

# --- 4.5: Ocena wydajności dla kryterium "co najmniej N narzędzi" ---
atLeastNPerformance <- sumAtLeast(sumSingleToolsList = singleToolPerformance)
bestToolAtLeastResults <- bestToolAtLeast(sumAtLeastList = atLeastNPerformance)
print("Wydajność dla kryterium 'co najmniej N narzędzi':")
print(bestToolAtLeastResults)

# --- 4.6: Ocena wydajności dla wszystkich kombinacji narzędzi ---
combPerformance <- sumCombTools(sumSingleToolsList = singleToolPerformance)
bestToolCombResults <- bestToolComb(sumCombToolsList = combPerformance)
print("Wydajność dla wybranych kombinacji narzędzi (top 5):")
print(head(bestToolCombResults, 5))

# --- 4.7: Oblicz macierze pomyłek i zwizualizuj wyniki ---
# Oblicz macierze pomyłek tylko dla kombinacji z precyzją > 0.9
allCms <- calculateCm(
  sumCombToolsList = combPerformance,
  bestToolCombMetricsData = bestToolCombResults, # Użyj wyników z bestToolComb
  printMetricThresholdMethods = TRUE,
  threshold = 0.9,
  returnOnlyHighMethods = TRUE,
  metricToExtract = "Precision"
)

# Wizualizacje (przekierowane do plików tymczasowych, aby uniknąć otwierania okien)
if (length(allCms) > 0) {
  png(tempfile()) # Przekieruj output graficzny
  radarPlotCm(cmList = allCms, layout = "multiple", displayArea = TRUE)
  clockPlotCm(cmList = allCms, layout = "multiple")
  dev.off() # Zamknij urządzenie graficzne
  message("Wygenerowano wykresy radarowe i zegarowe (do plików tymczasowych).")
}


# ===================================================================
# Sekcja 5: Analiza Funkcjonalna i Interakcji
# ===================================================================
message("Sekcja 5: Analiza funkcjonalna i interakcji...")

# --- 5.1: Zidentyfikuj ostateczną listę lncRNA ---
# Użyjemy wyników z `sumSingleTools`, aby znaleźć transkrypty przewidziane jako niekodujące
# przez co najmniej 5 narzędzi (wysoka pewność).
sumsDf <- as.data.frame(do.call(cbind, singleToolPerformance$tools))
predictedLncRnaIds <- singleToolPerformance$seqIDs[rowSums(sumsDf) >= 5]
# Połącz z listą znanych lncRNA
finalLncRnaTranscriptIds <- unique(c(predictedLncRnaIds, knownLncRnaIds))
finalLncRnaGeneIds <- unique(stringtieGtf[stringtieGtf$transcript_id %in% finalLncRnaTranscriptIds]$gene_id)
message(paste("Finalna liczba genów lncRNA do analizy:", length(finalLncRnaGeneIds)))

# --- 5.2: Analiza interakcji Cis ---
cisInteractionTable <- cisInter(
  feelncClasses = file.path(dataDir, "lncCodPot_MMus/FEELnc_Mm_classes.txt"),
  lncRNAs = finalLncRnaTranscriptIds,
  lncRnaLevel = "transcript",
  mRnaLevel = "gene",
  maxDist = 10000,
  mRNAs = rownames(degs)
)
if(nrow(cisInteractionTable) > 0) {
  cisGostResults <- gprofiler2::gost(query = unique(cisInteractionTable$targetRNAId), organism = "mmusculus")
  if(!is.null(cisGostResults)) {
    cisInteractionsProcessed <- processInteractions(gprof = cisGostResults$result, interaction_table = cisInteractionTable, type = "cis")
  }
}


# --- 5.3: Analiza interakcji Trans ---
transInteractionTable <- calculateTransInteractions(
  exprMatrix = genesCounts,
  rValueCutoff = 0.95,
  lncRnaIds = rownames(dels),
  targetRnaIds = rownames(degs)
)
if(nrow(transInteractionTable) > 0) {
  transGostResults <- gprofiler2::gost(query = unique(transInteractionTable$targetRNAId), organism = "mmusculus")
  if(!is.null(transGostResults)) {
    transInteractionsProcessed <- processInteractions(gprof = transGostResults$result, interaction_table = transInteractionTable, type = "trans")
  }
}


# --- 5.4: Połącz wszystkie przetworzone tabele interakcji ---
combinedInteractionsTable <- do.call(rbind, list(
  if(exists("cisInteractionsProcessed")) cisInteractionsProcessed,
  if(exists("transInteractionsProcessed")) transInteractionsProcessed
  # Tutaj można dodać wyniki z LncTar i LION, jeśli były obliczone
))


# --- 5.5: Wizualizacja wyników funkcjonalnych (Sankey Plots) ---
if (exists("combinedInteractionsTable") && nrow(combinedInteractionsTable) > 0) {
  # Przykładowe wizualizacje - odkomentuj, aby uruchomić
  # plotByTerms(data = combinedInteractionsTable, selectedTerms = "response to stress")
  # plotByType(data = combinedInteractionsTable, interactionTypes = "cis")
  message("Wykresy Sankeya można teraz generować na podstawie 'combinedInteractionsTable'.")
}

message("Analiza zakończona.")
# ===================================================================
# KONIEC SKRYPTU
# ===================================================================
