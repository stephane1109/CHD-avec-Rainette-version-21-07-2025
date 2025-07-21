################################
# # Rainette ## développé par Julien Barnier
# https://cran.r-project.org/web/packages/rainette/vignettes/introduction_usage.html
# https://juba.r-universe.dev/builds
################################


# install.packages(c("rainette", "quanteda", "wordcloud", "RColorBrewer", "dplyr", "shiny"))
# install.packages("htmltools")

###############################################################################
#                             Script CHD                                      #
#      A partir d'un corpus texte formaté aux exigences IRAMUTEQ              #
#                     version : 20-07-2025                                    #
#                                                                             #
#      1.Aide au paramétrage dans R de la CHD sur le corpus                   #
#      2.Extrait chi2, rainette_explor                                        #
#      3.Génère nuages de mots par classe                                     #
#      4.Exporte les segments de texte par classe au format text              #
#      5.Creation d'un concordancier au format html                           #
#      6.Affichage de la CHD avec rainette_explor (navigateur)                #
###############################################################################


# Chargement des bibliothèques nécessaires
library(rainette)
library(quanteda)
library(wordcloud)
library(RColorBrewer)
library(dplyr)
library(htmltools)
library(udpipe)

# Chemin du modèle UDPipe français
modele_udpipe <- "french-gsd-ud-2.5-191206.udpipe"
if (!file.exists(modele_udpipe)) {
  cat("Téléchargement du modèle UDPipe français...\n")
  udpipe_download_model(language = "french", model_dir = ".", overwrite = FALSE)
}
ud_model <- udpipe_load_model(modele_udpipe)

#########################################################
# PARAMÈTRES UTILISATEUR (modifiable par l'utilisateur)
#########################################################

# Taille des segments avant analyse (nombre de mots par segment)
segment_size <- 40

# Nom du fichier texte brut (dans base_dir)
texte_fichier <- "tourisme_surtourisme_1an_presse_nationale.txt"

# Répertoire principal contenant le fichier texte
base_dir <- "/Users/stephanemeurisse/Documents/Cours ISTHIA/Cours ISTHIA 2026/Cours_Analyse_textuelle/"

# Nombre de classes (clusters) pour Rainette
k <- 4

# Nombre minimal de segments de texte par classe pour continuer à diviser
# (évite de créer des classes trop petites)
min_split_segments <- 10

# Seuil minimal de fréquence documentaire dans le DFM
# Exemple : si min_docfreq <- 2 seuls les mots présents dans au moins 2 segments sont gardés
# Les mots présents dans seulement 1 segment sont retirés
min_docfreq <- 2

# Nombre maximal de mots affichés par classe dans les nuages de mots
top_n <- 20

# Répertoire où exporter les résultats
export_dir <- file.path(base_dir, "exports")
dir.create(export_dir, showWarnings = FALSE, recursive = TRUE)

#########################################################
# Chargement et découpage du corpus
#########################################################
chemin_fichier <- file.path(base_dir, texte_fichier)
corpus <- import_corpus_iramuteq(chemin_fichier)
cat("Nombre de documents importés :", ndoc(corpus), "\n")

corpus <- split_segments(corpus, segment_size = segment_size)
cat("Nombre de segments après découpage :", ndoc(corpus), "\n")

docvars(corpus, "orig_doc_id") <- gsub("_.*", "", docnames(corpus))

#########################################################
# Lemmatisation avec UDPipe
#########################################################
cat("Lemmatisation en cours...\n")
annotation_list <- lapply(seq_along(corpus), function(i) {
  cat("  Segment", i, "/", ndoc(corpus), "\n")
  res <- udpipe_annotate(ud_model, x = as.character(corpus[i]))
  df <- as.data.frame(res)
  df$doc_id <- docnames(corpus)[i]
  return(df)
})
annotation_df <- do.call(rbind, annotation_list)

annotation_df <- annotation_df %>%
  filter(upos %in% c("NOUN", "ADJ", "VERB"), !is.na(lemma), lemma != "", lemma != " ")

cat("Nombre total de lemmes significatifs :", nrow(annotation_df), "\n")

textes_lemmat <- annotation_df %>%
  group_by(doc_id) %>%
  summarise(text = paste(lemma, collapse = " ")) %>%
  ungroup()

cat("Nombre unique de doc_id après regroupement :", length(unique(textes_lemmat$doc_id)), "\n")

corpus <- corpus(textes_lemmat$text, docnames = textes_lemmat$doc_id)
docvars(corpus, "orig_doc_id") <- gsub("_.*", "", docnames(corpus))

cat("Extraits des textes lemmatisés :\n")
for (i in 1:min(5, ndoc(corpus))) {
  cat(paste0(">> ", docnames(corpus)[i], " :\n"))
  cat(substr(as.character(corpus)[i], 1, 300), "\n\n")
}

#########################################################
# Prétraitement du texte et création du DFM
#########################################################
tok <- tokens(corpus, remove_punct = TRUE, remove_numbers = TRUE)
tok <- tokens_split(tok, "'")
tok <- tokens_remove(tok, stopwords("fr"))
tok <- tokens_tolower(tok)

dfm <- dfm(tok)
dfm <- dfm_trim(dfm, min_docfreq = min_docfreq)

docvars(dfm, "doc_id") <- docvars(corpus, "orig_doc_id")

cat("Dimensions du DFM après trim :", dim(dfm), "\n")

included_segments <- docnames(dfm)
filtered_corpus <- corpus[included_segments]
docvars(filtered_corpus, "orig_doc_id") <- docvars(corpus[included_segments], "orig_doc_id")

cat("Vérification finale des doc_id avant rainette :\n")
print(head(docnames(dfm)))
print(head(docvars(dfm, "doc_id")))

#########################################################
# Classification hiérarchique descendante (CHD)
#########################################################
res <- rainette(dfm, k = k, min_segment_size = 0, min_split_members = min_split_members)
docvars(filtered_corpus)$Classes <- res$group

cat("DEBUG : Table de répartition des classes :\n")
print(table(docvars(filtered_corpus)$Classes))

#########################################################
# Extraction et export des segments par classe
#########################################################
segments_by_class <- split(as.character(filtered_corpus), docvars(filtered_corpus)$Classes)

segments_file <- file.path(export_dir, "segments_par_classe.txt")
writeLines(
  unlist(lapply(
    names(segments_by_class),
    function(cl) {
      c(paste0("Classe ", cl, ":"), segments_by_class[[cl]], "")
    }
  )),
  segments_file
)
cat("Segments par classe exportés dans :", segments_file, "\n")

#########################################################
# Récupération des termes discriminants (chi²)
#########################################################
res_stats_list <- rainette_stats(
  dtm = dfm,
  groups = docvars(filtered_corpus)$Classes,
  measure = c("chi2"),
  n_terms = 9999,
  show_negative = TRUE,
  max_p = 0.05
)

res_stats_df <- bind_rows(res_stats_list, .id = "Classe")

highlight_terms <- res_stats_df %>%
  filter(p <= 0.05, nchar(feature) >= 3) %>%
  group_by(Classe) %>%
  arrange(desc(chi2)) %>%
  slice_head(n = top_n) %>%
  pull(feature) %>%
  unique()

#########################################################
# Export HTML avec surlignage
#########################################################
html_file <- file.path(export_dir, "segments_par_classe.html")

highlight_text_html <- function(text, terms, start_tag, end_tag) {
  for (term in terms) {
    escaped_term <- gsub("([\\^\\$\\*\\+\\?\\(\\)\\[\\]\\{\\}\\.\\|])", "\\\\\\1", term)
    text <- gsub(
      paste0("\\b", escaped_term, "\\b"),
      paste0(start_tag, term, end_tag),
      text,
      ignore.case = TRUE
    )
  }
  return(text)
}

if (file.exists(html_file)) file.remove(html_file)

cat("<html><head><style>body { font-family: Arial; } span.highlight { background-color: yellow; }</style></head><body>\n",
    file = html_file, append = TRUE)
cat("<h1>Segments par Classe (avec surlignage des termes discriminants)</h1>\n", file = html_file, append = TRUE)

for (cl in names(segments_by_class)) {
  cat(paste0("<h2>Classe ", cl, "</h2>\n"), file = html_file, append = TRUE)
  for (segment in segments_by_class[[cl]]) {
    highlighted_segment <- tryCatch(
      highlight_text_html(segment, highlight_terms, "<span class='highlight'>", "</span>"),
      error = function(e) segment
    )
    cat(paste0("<p>", highlighted_segment, "</p>\n"), file = html_file, append = TRUE)
  }
}

cat("</body></html>\n", file = html_file, append = TRUE)
cat("Fichier HTML exporté dans :", html_file, "\n")

#########################################################
# Visualisations : Nuages de mots
#########################################################
wordcloud_dir <- file.path(export_dir, "wordclouds")
dir.create(wordcloud_dir, showWarnings = FALSE, recursive = TRUE)

clusters <- sort(unique(docvars(filtered_corpus)$Classes))
for (cl in clusters) {
  subset_stats <- subset(res_stats_df, Classe == cl & p <= 0.05)
  subset_stats <- subset_stats[order(-subset_stats$chi2), ]
  
  if (nrow(subset_stats) > 0) {
    subset_stats <- head(subset_stats, top_n)
    
    png_filename <- file.path(wordcloud_dir, paste0("cluster_", cl, "_wordcloud.png"))
    png(png_filename, width = 800, height = 600)
    
    wordcloud(
      words = subset_stats$feature,
      freq  = subset_stats$chi2,
      scale = c(10, 0.5),
      max.words = top_n,
      colors = brewer.pal(8, "Dark2")
    )
    dev.off()
    cat("Nuage de mots pour la classe", cl, "enregistré dans :", png_filename, "\n")
  } else {
    cat("Classe", cl, ": aucun terme significatif pour le nuage de mots.\n")
  }
}

#########################################################
# Affichage interactif rainette
#########################################################
rainette_explor(res, dfm, filtered_corpus)
