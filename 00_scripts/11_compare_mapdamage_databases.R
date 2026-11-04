#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE, warn = 1)

# ============================================================
# Configuration
# ============================================================

args <- commandArgs(trailingOnly = TRUE)

base_dir <- if (length(args) >= 1L) {
  args[1L]
} else {
  paste0(
    "/home/plstenge/seda_DNA_Corsican_wreck_comparison/",
    "10_mapdamage_cascade/runs/run_20261103_150905_9459089/results"
  )
}

# Sorties hors du dossier results, pour ne pas les confondre
# avec les dossiers d'especes lors d'une execution ulterieure.
out_dir <- if (length(args) >= 2L) {
  args[2L]
} else {
  file.path(dirname(base_dir), "comparison_mapdamage_databases")
}

samples <- c("sed6", "sed8")
databases <- c("k25", "k29", "k35", "k35_k29_k25")
highlight_db <- "k35_k29_k25"

# Nombre de positions affichees depuis chaque extremite.
n_positions <- 25L

# TRUE : meme echelle Y pour sed6 et sed8 d'une meme espece.
# FALSE : echelle propre a chaque graphique.
same_y_per_species <- TRUE

# Frequences exprimees entre 0 et 1, comme sur la figure fournie.
# Aucun lissage ni moyenne entre bases ou echantillons.
db_colors <- c(
  k25         = "#E8BD8F",
  k29         = "#A9BBD5",
  k35         = "#B5CDB1",
  k35_k29_k25 = "#B2182B"
)

db_lwd <- c(
  k25         = 1.3,
  k29         = 1.3,
  k35         = 1.3,
  k35_k29_k25 = 3.2
)

db_cex <- c(
  k25         = 0.40,
  k29         = 0.40,
  k35         = 0.40,
  k35_k29_k25 = 0.65
)

# ============================================================
# Verification des chemins
# ============================================================

if (!dir.exists(base_dir)) {
  stop("Dossier de resultats introuvable : ", base_dir)
}

base_dir <- normalizePath(base_dir, mustWork = TRUE)

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

if (!dir.exists(out_dir)) {
  stop("Impossible de creer le dossier de sortie : ", out_dir)
}

species_dirs <- list.dirs(
  base_dir,
  full.names = TRUE,
  recursive = FALSE
)

# Ne retenir que les dossiers contenant sed6 et/ou sed8.
species_dirs <- species_dirs[
  vapply(
    species_dirs,
    function(d) any(dir.exists(file.path(d, samples))),
    logical(1)
  )
]

species_dirs <- sort(species_dirs)

if (length(species_dirs) == 0L) {
  stop("Aucun dossier d'espece contenant sed6 ou sed8.")
}

# ============================================================
# Lecture des fichiers de frequences
# ============================================================

empty_profile <- function() {
  data.frame(
    position = integer(),
    frequency = numeric()
  )
}

read_damage <- function(path, max_position) {

  if (!file.exists(path)) {
    warning("Fichier absent : ", path, call. = FALSE)
    return(empty_profile())
  }

  txt <- tryCatch(
    readLines(path, warn = FALSE),
    error = function(e) {
      warning(
        "Lecture impossible : ", path, " : ", conditionMessage(e),
        call. = FALSE
      )
      character()
    }
  )

  # Supprimer lignes vides et commentaires.
  txt <- trimws(txt)
  txt <- txt[nzchar(txt) & !grepl("^#", txt)]

  if (length(txt) == 0L) {
    warning("Fichier vide : ", path, call. = FALSE)
    return(empty_profile())
  }

  # Lecture des deux premieres colonnes, separees par espaces
  # ou tabulations. Les en-tetes non numeriques sont ignores.
  fields <- strsplit(txt, "[[:space:]]+")

  position <- vapply(
    fields,
    function(z) suppressWarnings(as.numeric(z[1L])),
    numeric(1)
  )

  frequency <- vapply(
    fields,
    function(z) {
      if (length(z) < 2L) return(NA_real_)
      suppressWarnings(as.numeric(z[2L]))
    },
    numeric(1)
  )

  valid_position <- is.finite(position) &
    position >= 1 &
    position == floor(position)

  if (!any(valid_position)) {
    warning(
      "Aucune position positive lisible dans : ", path,
      call. = FALSE
    )
    return(empty_profile())
  }

  # Ne pas transformer les valeurs manquantes en zeros.
  invalid_frequency <- !is.finite(frequency) |
    frequency < 0 | frequency > 1

  n_invalid <- sum(valid_position & invalid_frequency)

  if (n_invalid > 0L) {
    warning(
      n_invalid, " frequence(s) absente(s) ou invalide(s) dans : ",
      path, " ; conservees comme NA.",
      call. = FALSE
    )
  }

  frequency[invalid_frequency] <- NA_real_

  keep <- valid_position & position <= max_position

  d <- data.frame(
    position = as.integer(position[keep]),
    frequency = frequency[keep]
  )

  if (anyDuplicated(d$position)) {
    stop(
      "Positions dupliquees dans : ", path,
      "\nVerifier le format du fichier avant de tracer."
    )
  }

  d <- d[order(d$position), , drop = FALSE]

  if (!any(is.finite(d$frequency))) {
    warning(
      "Aucune frequence exploitable entre les positions 1 et ",
      max_position, " : ", path,
      call. = FALSE
    )
  }

  d
}

# ============================================================
# Chargement de toutes les donnees
# ============================================================

rows <- list()

for (species_dir in species_dirs) {

  species_id <- basename(species_dir)

  for (sample_id in samples) {

    sample_dir <- file.path(species_dir, sample_id)

    if (!dir.exists(sample_dir)) {
      warning("Echantillon absent : ", sample_dir, call. = FALSE)
      next
    }

    for (db in databases) {

      md_dir <- file.path(sample_dir, db, "mapDamage")

      for (end_id in c("5p", "3p")) {

        filename <- if (end_id == "5p") {
          "5pCtoT_freq.txt"
        } else {
          "3pGtoA_freq.txt"
        }

        d <- read_damage(
          file.path(md_dir, filename),
          max_position = n_positions
        )

        if (nrow(d) == 0L) next

        rows[[length(rows) + 1L]] <- data.frame(
          species = species_id,
          sample = sample_id,
          database = db,
          end = end_id,
          substitution = if (end_id == "5p") "C>T" else "G>A",
          position = d$position,
          frequency = d$frequency
        )
      }
    }
  }
}

if (length(rows) == 0L) {
  stop("Aucun profil lu. Verifier les chemins et fichiers d'entree.")
}

damage <- do.call(rbind, rows)

if (!any(is.finite(damage$frequency))) {
  stop("Aucune frequence exploitable dans les fichiers lus.")
}

write.table(
  damage,
  file = file.path(out_dir, "mapdamage_profiles_1_to_25.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  na = "NA"
)

# ============================================================
# Fonctions graphiques
# ============================================================

get_ymax <- function(d) {
  y <- d$frequency[is.finite(d$frequency)]
  if (length(y) == 0L) return(0.05)
  min(1, max(0.05, max(y) * 1.15))
}

plot_profile <- function(d, title, ymax) {

  # 5' : positions 1 -> 25.
  # 3' : positions 25 -> 1.
  # Un espace separe les deux moities.
  gap <- 2L
  right_start <- n_positions + gap + 1L
  right_end <- 2L * n_positions + gap
  separator <- (n_positions + right_start) / 2

  par(
    mar = c(5.1, 5.0, 4.0, 10.8),
    mgp = c(2.9, 0.8, 0),
    las = 1,
    family = "sans"
  )

  plot(
    NA,
    xlim = c(0.5, right_end + 0.5),
    ylim = c(0, ymax),
    xaxs = "i",
    yaxs = "i",
    xaxt = "n",
    yaxt = "n",
    bty = "n",
    xlab = "Position depuis les extremites du fragment (nt)",
    ylab = "Frequence de misincorporation"
  )

  yticks <- pretty(c(0, ymax), n = 5)
  yticks <- yticks[yticks >= 0 & yticks <= ymax]

  abline(h = yticks, col = "#EEEEEE", lwd = 0.8)
  abline(v = separator, col = "#BBBBBB", lty = 3)

  axis(2, at = yticks, cex.axis = 0.9)
  axis(
    1,
    at = c(seq_len(n_positions), right_start:right_end),
    labels = c(seq_len(n_positions), n_positions:1L),
    cex.axis = 0.66,
    tck = -0.018
  )

  # Les courbes secondaires sont tracees en premier ;
  # la combinaison est donc toujours tracee par-dessus.
  draw_order <- c(setdiff(databases, highlight_db), highlight_db)

  for (db in draw_order) {
    for (end_id in c("5p", "3p")) {

      z <- d[
        d$database == db & d$end == end_id,
        c("position", "frequency"),
        drop = FALSE
      ]

      # Grille complete pour interrompre la ligne si une
      # position est absente, plutot que de relier au travers.
      y <- rep(NA_real_, n_positions)

      if (nrow(z) > 0L) {
        y[z$position] <- z$frequency
      }

      x <- if (end_id == "5p") {
        seq_len(n_positions)
      } else {
        right_end + 1L - seq_len(n_positions)
      }

      lines(
        x, y,
        col = db_colors[[db]],
        lwd = db_lwd[[db]]
      )

      points(
        x, y,
        pch = 16,
        col = db_colors[[db]],
        cex = db_cex[[db]]
      )
    }
  }

  text(
    1, ymax * 0.95, "5' C>T",
    adj = c(0, 1), font = 2, cex = 1
  )

  text(
    right_end, ymax * 0.95, "3' G>A",
    adj = c(1, 1), font = 2, cex = 1
  )

  title(main = title, adj = 0, cex.main = 1.05)

  mtext(
    paste0(
      "Positions 1-", n_positions,
      " ; combinaison mise en evidence"
    ),
    side = 3, line = 0.35, adj = 0, cex = 0.8
  )

  available_db <- databases[
    vapply(
      databases,
      function(db) {
        any(d$database == db & is.finite(d$frequency))
      },
      logical(1)
    )
  ]

  legend(
    x = right_end + 2,
    y = ymax * 0.78,
    legend = available_db,
    title = "Database",
    col = db_colors[available_db],
    lwd = db_lwd[available_db],
    lty = 1,
    pch = 16,
    pt.cex = db_cex[available_db],
    bty = "n",
    xpd = NA,
    cex = 0.88,
    y.intersp = 1.3
  )

  # Signaler explicitement les profils absents ou illisibles.
  missing_profiles <- character()

  for (db in databases) {
    for (end_id in c("5p", "3p")) {
      ok <- any(
        d$database == db &
          d$end == end_id &
          is.finite(d$frequency)
      )

      if (!ok) {
        missing_profiles <- c(
          missing_profiles,
          paste(db, end_id, sep = "/")
        )
      }
    }
  }

  if (length(missing_profiles) > 0L) {
    note <- paste(
      strwrap(
        paste(
          "Profils absents :",
          paste(missing_profiles, collapse = ", ")
        ),
        width = 105
      ),
      collapse = "\n"
    )

    mtext(
      note,
      side = 1,
      line = 4.0,
      cex = 0.62,
      col = "#A00000"
    )
  }
}

save_pdf <- function(path, d, title, ymax) {
  pdf(
    path,
    width = 11,
    height = 6,
    useDingbats = FALSE
  )
  on.exit(dev.off(), add = TRUE)
  plot_profile(d, title, ymax)
}

save_png <- function(path, d, title, ymax) {

  # Cairo permet l'export PNG sans serveur X sur le cluster.
  if (!capabilities("cairo")) {
    warning(
      "Support Cairo indisponible : PNG ignore, PDF conserve.",
      call. = FALSE
    )
    return(invisible(NULL))
  }

  png(
    path,
    width = 3300,
    height = 1800,
    res = 300,
    type = "cairo"
  )
  on.exit(dev.off(), add = TRUE)
  plot_profile(d, title, ymax)
}

# ============================================================
# Une figure par espece et par echantillon
# ============================================================

jobs <- list()

for (species_dir in species_dirs) {

  species_id <- basename(species_dir)

  # Etiquette courte, en conservant l'identifiant complet
  # de l'assemblage sur une deuxieme ligne.
  parts <- strsplit(species_id, "_", fixed = TRUE)[[1L]]

  species_label <- if (length(parts) >= 2L) {
    paste(parts[1:2], collapse = " ")
  } else {
    species_id
  }

  species_data <- damage[
    damage$species == species_id,
    ,
    drop = FALSE
  ]

  species_ymax <- get_ymax(species_data)

  species_out <- file.path(out_dir, species_id)
  dir.create(species_out, recursive = TRUE, showWarnings = FALSE)

  for (sample_id in samples) {

    d <- species_data[
      species_data$sample == sample_id,
      ,
      drop = FALSE
    ]

    if (!any(is.finite(d$frequency))) {
      warning(
        "Pas de figure pour ", species_id, " / ", sample_id,
        " : aucune frequence exploitable.",
        call. = FALSE
      )
      next
    }

    ymax <- if (same_y_per_species) {
      species_ymax
    } else {
      get_ymax(d)
    }

    title <- paste0(
      sample_id, " - ", species_label,
      "\n", species_id
    )

    prefix <- file.path(
      species_out,
      paste0(species_id, "__", sample_id, "__mapdamage_databases")
    )

    save_pdf(paste0(prefix, ".pdf"), d, title, ymax)
    save_png(paste0(prefix, ".png"), d, title, ymax)

    jobs[[length(jobs) + 1L]] <- list(
      data = d,
      title = title,
      ymax = ymax
    )

    message("Figure creee : ", species_id, " / ", sample_id)
  }
}

# ============================================================
# PDF multipage : sed6 puis sed8 pour chaque espece
# ============================================================

save_combined_pdf <- function(path, jobs) {
  pdf(
    path,
    width = 11,
    height = 6,
    onefile = TRUE,
    useDingbats = FALSE
  )
  on.exit(dev.off(), add = TRUE)

  for (job in jobs) {
    plot_profile(job$data, job$title, job$ymax)
  }
}

if (length(jobs) > 0L) {
  save_combined_pdf(
    file.path(out_dir, "all_species_sed6_sed8_mapdamage.pdf"),
    jobs
  )
}

cat(
  "\nTermine : ", length(jobs), " graphiques.\n",
  "Dossier de sortie : ", out_dir, "\n",
  sep = ""
)