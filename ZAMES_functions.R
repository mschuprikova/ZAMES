
# ==============================================================================
# -----------------------------------Packages-----------------------------------
# ==============================================================================

library("sangerseqR")   # import of ABIF files
library("Biostrings")   # import of FASTA and character inputs
library("stringr")      # processing of channel order (and filepaths)

# ==============================================================================
# -------------------------------Data Import Step-------------------------------
# ==============================================================================

# ----------------- Function 1: import of date and time info -------------------
dateTime <- function(
    date, 
    time
    ) {
  
  result <- as.POSIXct(
    sprintf(
      "%04d-%02d-%02d %02d:%02d:%02d.%02d",
      date$year,
      date$month,
      date$day,
      time$hour,
      time$minute,
      time$second,
      time$hsecond
    ),
    format = "%Y-%m-%d %H:%M:%OS",
    tz = ""
  )
  
  return(result)
}

# ----------------- Function 2: import of date and time info -------------------
dateTimeImport <- function(
    date_start, 
    date_stop, 
    time_start, 
    time_stop) {
  
  start <- dateTime(
    date = date_start, 
    time = time_start
    )
  
  stop <- dateTime(
    date = date_stop, 
    time = time_stop
    )
  
  duration <- as.POSIXct(stop, format = "%Y-%m-%d %H:%M:%OS") - as.POSIXct(start, format = "%Y-%m-%d %H:%M:%OS")
  duration_s <- round(as.numeric(duration, units = "secs"))
  duration_hms <- sprintf(
    "%02d:%02d:%02d",
    duration_s %/% 3600,
    (duration_s %% 3600) %/% 60,
    duration_s %% 60
  )
  
  info <- data.frame(
    start_date = as.Date(start),
    start_time = format(start, "%H:%M:%OS"),
    start = start,
    
    stop_date = as.Date(stop),
    stop_time = format(stop, "%H:%M:%OS"),
    stop = stop,
    
    duration = duration_hms
  )
  
  return(info)
}

# ------------------ Function 3: import of one sanger file ---------------------
readImport <- function(read, filename) {
  
  require("sangerseqR")
  
  # Импорт и подготовка
  raw <- read.abif(read)   # raw ABIF data
  raw_obj <- sangerseq(raw)   # Need for obtain more decoded seq

  # Получение основной информации
  seq <- primarySeq(raw_obj)   # DNAString Object
  basecalled <- makeBaseCalls(raw_obj)   # DNAString Object
  seq_basecalled <- primarySeq(basecalled)
  pos <- raw@data$PLOC.1    # list of numeric
  qual <- raw@data$PCON.1   # list of numeric
  sigs <- traceMatrix(raw_obj)
  
  # Получение дополнительной информации
  run_date_start <- tryCatch(
    raw@data$RUND.1, 
    error = function(e) NA
    )
  run_date_stop <- tryCatch(
    raw@data$RUND.2, 
    error = function(e) NA
    )
  run_time_start <- tryCatch(
    raw@data$RUNT.1, 
    error = function(e) NA
    )
  run_time_stop <- tryCatch(
    raw@data$RUNT.2, 
    error = function(e) NA
    )
  
  scan_date_start <- tryCatch(
    raw@data$RUND.3, 
    error = function(e) NA
    )
  scan_date_stop <- tryCatch(
    raw@data$RUND.4, 
    error = function(e) NA
    )
  scan_time_start <- tryCatch(
    raw@data$RUNT.3, 
    error = function(e) NA
    )
  scan_time_stop <- tryCatch(
    raw@data$RUNT.4, 
    error = function(e) NA
    )
  
  run_info <- dateTimeImport(
    date_start = run_date_start, 
    date_stop = run_date_stop, 
    time_start = run_time_start, 
    time_stop = run_time_stop
    )
  
  scan_info <- dateTimeImport(
    date_start = scan_date_start, 
    date_stop = scan_date_stop, 
    time_start = scan_time_start, 
    time_stop = scan_time_stop
  )
  
  user <- tryCatch(
    raw@data$User.1, 
    error = function(e) NA
    )
  machine <- tryCatch(
    raw@data$MCHN, 
    error = function(e) NA
    )
  dye_kit <- tryCatch(
    raw@data$DySN.1, 
    error = function(e) NA
    )
    
  # Вывод
  main_info <- list(
    filename = filename,
    seq = seq,
    seq_basecalled = seq_basecalled,
    pos = pos,
    qual = qual,
    sigs = sigs
  )
  
  additional_info <- list(
    run_info = run_info,
    scan_info = scan_info,
    user = user,
    machine = machine,
    dye_kit = dye_kit
  )
  
  result <- list(
    main_info = main_info,
    additional_info = additional_info
  )
  
  return(result)
}

# ----------- Function 4: run import of one sanger file for batch --------------
readImportRun <- function(
    batch, 
    filenames
    ) {
  
  imported <- list()
  
  for (i in 1:length(batch)) {
    
    read_info <- readImport(
      read = batch[i],
      filename = filenames[[i]]
    )
    
    imported[[filenames[[i]]]] <- read_info
  }

  return(imported)
}

# ------------- Function 5: import of target and reference info ----------------
referenceTargetImport <- function(
    reference, 
    target_region, 
    target_site
    ){
  
  require("Biostrings")
  
  ref_seq <- readDNAStringSet(reference)
  target_seq <- DNAString(target_region)
  
  target_info <- list(
    target_seq = target_seq,
    target_site = target_site
  )
  result <- list(
    reference_info = ref_seq,
    target_info = target_info
  )
  return(result)
}

# -------------------- Function 6: data import main run ------------------------
dataImportRun <- function(input) {
  
  control_imported_info <- readImportRun(
    batch = input$control_files$datapath,
    filenames = input$control_files$name
    )
  
  sample_imported_info <- readImportRun(
    batch = input$sample_files$datapath,
    filenames = input$sample_files$name
    )
  
  reference_target_imported_info <- referenceTargetImport(
    reference = input$reference_file$datapath,
    target_region = input$target_region,
    target_site = input$target_site
  )
  
  imported <- list(
    control_info = control_imported_info,
    sample_info = sample_imported_info,
    reference_target_info = reference_target_imported_info
  )
  
  return(imported)
}

# ==============================================================================
# ---------------------------Direction identifying Step-------------------------
# ==============================================================================

# ----------------------- Function 1: Score identifying ------------------------
scoreIdentifying <- function(
    pattern, 
    subject
    ) {
  
  require("Biostrings")
  
  # Выравнивание subject на pattern
  aln_fwd <- pairwiseAlignment(
    pattern = pattern, 
    subject = subject, 
    type = "local"
    )
  
  aln_rev <- pairwiseAlignment(
    pattern = reverseComplement(pattern), 
    subject = subject, 
    type = "local"
    )
  
  # Проверка на совпадение
  fwd_identical <- identical(
    as.character(pattern), 
    as.character(subject(aln_fwd))
    )
  
  rev_identical <- identical(
    as.character(reverseComplement(pattern)), 
    as.character(subject(aln_rev))
    )
  
  if (!fwd_identical && !rev_identical) {
    warning("⚠️  Subject does not perfectly match reference in either orientation. Please verify pattern or subject sequence.")
  }
  if (fwd_identical && rev_identical) {
    warning("⚠️  Subject matches the pattern in both directions. Please verify subject sequence.")
  }
  
  score_fwd <- score(aln_fwd)
  score_rev <- score(aln_rev)
  
  result <- list(
    score_fwd = score_fwd,
    score_rev = score_rev,
    aln_fwd = aln_fwd,
    aln_rev = aln_rev
  )
  
  return(result)
}

# ---------------------- Function 2: Sequence direction ------------------------
sequenceDirection <- function(score_fwd, score_rev, aln_fwd, aln_rev,
                               mode, ref_score = NULL, identity_low = 0.1) {
  
  if (mode == "ref") {
    if (score_fwd >= score_rev) {
      direction <- "Forward"
      aln <- aln_fwd
      score <- score_fwd
    } else if (score_fwd <= score_rev) {
      direction <- "Reverse"
      aln <- aln_rev
      score <- score_rev
    } else {
      direction <- "Undetermined"
      warning("⚠️  The match is weak or ambiguous. Possibly due to mutations or sequencing noise.")
    }
    
    result <- list(
      direction = direction,
      ref_score = score
    )
  }
  
  if (mode == "read") {
    if (score_fwd >= score_rev && score_fwd > ref_score * identity_low) {
      direction <- "Forward"
      aln <- aln_fwd
    } else if (score_rev > score_fwd && score_rev > ref_score * identity_low) {
      direction <- "Reverse"
      aln <- aln_rev
    } else {
      direction <- "Undetermined"
      aln <- if (score_fwd >= score_rev) aln_fwd else aln_rev
      warning("⚠️  The match is weak or ambiguous. Possibly due to mutations or sequencing noise.")
    }
    result <- list(
      direction = direction,
      aln = aln
    )
  }
  
  return(result)
}

# ----------------------- Function 3: Target boundaries ------------------------
targetBoundaries <- function(aln, direction, ref_score, target_region, target_site, identity_high = 0.9) {
  
  require("Biostrings")
  
  # Полное совпадение
  if (score(aln) == ref_score) {
    start_pos <- start(subject(aln))
    end_pos <- end(subject(aln))
    if (direction == "Forward") {
      ts_pos <- start_pos + target_site - 1
      }
    else if (direction == "Reverse") {
      ts_pos <- end_pos - target_site + 1
      }
    }

  # Неполное совпадение
  if (score(aln) < ref_score * identity_high) {
    # Границы Target Region
    if (start(pattern(aln)) == 1) {
      start_pos <- start(subject(aln))
      end_pos <- start_pos + str_length(as.character(target_region)) - 1
    } else if (start(pattern(aln)) >= 1) {
      end_pos <- end(subject(aln))
      start_pos <- end_pos - str_length(as.character(target_region)) + 1
    }
    
    # Позиция Target Site в Read
    if (direction == "Forward") {
      ts_pos <- start_pos + target_site - 1
    } else if (direction == "Reverse") {
      ts_pos <- end_pos - target_site + 1
      }
  }
  
  result <- list(
    direction = direction,
    start = start_pos,
    end = end_pos,
    ts = ts_pos
  )
  
  return(result)
}

# ---------------------------- Function 4: Read Run -----------------=----------
readDirectionRun <- function(reference, read, target_region, target_site) {
  
  require("Biostrings")

  ### --- 1 --- Поиск ref score
  ref_score_info <- scoreIdentifying(target_region, reference)
  
  ref_direction <- sequenceDirection(score_fwd = ref_score_info$score_fwd,
                                     score_rev = ref_score_info$score_rev,
                                     aln_fwd = ref_score_info$aln_fwd,
                                     aln_rev = ref_score_info$aln_rev,
                                     mode = "ref")
  
  ### --- 2 --- Reads direction
  read_score_info <- scoreIdentifying(target_region, read)
  
  read_direction_info <- sequenceDirection(score_fwd = read_score_info$score_fwd,
                                           score_rev = read_score_info$score_rev,
                                           aln_fwd = read_score_info$aln_fwd,
                                           aln_rev = read_score_info$aln_rev,
                                           ref_score = ref_direction$ref_score,
                                           mode = "read")
  
  target_info <- targetBoundaries(aln = read_direction_info$aln, 
                                  direction = read_direction_info$direction, 
                                  ref_score = ref_direction$ref_score, 
                                  target_region = target_region, 
                                  target_site = target_site)
  
  result <- list(
    target_info = target_info,
    ref_direction = ref_direction$direction
  )
  
  return(result)
}

# ---------------------------- Function 4: Batch Run -----------------=---------
directionRun <- function(imported) {
  
  control_batch <- imported$control_info
  sample_batch <- imported$sample_info
  reference <- imported$reference_target_info$reference_info
  target_region <- imported$reference_target_info$target_info$target_seq
  target_site <- imported$reference_target_info$target_info$target_site
  
  # Control
  for (i in seq_along(control_batch)) {
    read_direction_info <- readDirectionRun(
      reference = reference,
      read = control_batch[[i]]$main_info$seq,
      target_region = target_region,
      target_site = target_site
    )
    control_batch[[i]]$aln_info <- read_direction_info$target_info
    control_batch[[i]]$target_info$direction <- read_direction_info$ref_direction
  }
  
  # Sample
  for (i in seq_along(sample_batch)) {
    read_direction_info <- readDirectionRun(
      reference = reference,
      read = sample_batch[[i]]$main_info$seq,
      target_region = target_region,
      target_site = target_site
    )
    sample_batch[[i]]$aln_info <- read_direction_info$target_info
    sample_batch[[i]]$target_info$direction <- read_direction_info$ref_direction
  }
  
  # Возвращаем обновлённые batch в imported
  imported$control_info <- control_batch
  imported$sample_info <- sample_batch
  
  return(imported)
}

# ==============================================================================
# ---------------------------- Read conversion Step 3 --------------------------
# ==============================================================================

# -------------------------- Function 1: Reverse Read --------------------------
reverseRead <- function(read) {
  
  # Если уже Forward — ничего не делаем
  if (read$aln_info$direction == read$target_info$direction) {
    return(read)
  }
  
  # 1. Длина нуклеотидной последовательности
  seq_len <- length(read$main_info$seq)
  
  # 2. Reverse complement последовательностей
  read$main_info$seq <- reverseComplement(read$main_info$seq)
  read$main_info$seq_basecalled <- reverseComplement(read$main_info$seq_basecalled)
  
  # 3. Reverse качества
  read$main_info$qual <- rev(read$main_info$qual)
  
  # 4. Reverse позиции пиков
  trace_len <- nrow(read$main_info$sigs)
  read$main_info$pos <- trace_len + 1 - rev(read$main_info$pos)
  
  # 5. Reverse chromatogram
  sigs <- read$main_info$sigs
  sigs <- sigs[nrow(sigs):1, , drop = FALSE]   # reverse
  sigs <- sigs[, c(4, 3, 2, 1), drop = FALSE]  # complement
  read$main_info$sigs <- sigs
  
  # 6. Пересчёт координат alignment
  old_start <- read$aln_info$start
  old_end   <- read$aln_info$end
  old_ts    <- read$aln_info$ts
  
  read$aln_info$start <- seq_len + 1 - old_end
  read$aln_info$end   <- seq_len + 1 - old_start
  read$aln_info$ts    <- seq_len + 1 - old_ts
  
  return(read)
}

# ------------------ Function 2: Preparation to visualization ------------------
prepareRead <- function(read) {
  
  sigs <- read$main_info$sigs
  pos <- read$main_info$pos
  seq <- strsplit(as.character(read$main_info$seq), "")[[1]]
  filename <- read$main_info$filename
  
  # Проверки
  if (length(pos) != length(seq)) {
    stop("Длина pos не совпадает с длиной seq.")
  }
  
  if (ncol(sigs) != 4) {
    stop("sigs должен содержать 4 столбца: A, C, G, T.")
  }
  
  # Данные сигналов
  signal_df <- data.frame(
    trace_pos = seq_len(nrow(sigs)),
    A = sigs[, 1],
    C = sigs[, 2],
    G = sigs[, 3],
    T = sigs[, 4]
  )
  
  signal_df <- reshape(
    signal_df,
    varying = c("A", "C", "G", "T"),
    v.names = "signal",
    timevar = "nuc",
    times = c("A", "C", "G", "T"),
    direction = "long"
  )
  
  rownames(signal_df) <- NULL
  
  # Информация о пиках
  peak_df <- data.frame(nuc = seq,
                        peak_index = seq_along(seq),
                        peak_pos = pos)
  
  # aln_info
  start_index <- read$aln_info$start
  end_index   <- read$aln_info$end
  ts_index    <- read$aln_info$ts
  
  start_pos <- pos[start_index]
  end_pos   <- pos[end_index]
  ts_pos    <- pos[ts_index]
  
  result <- list(
    signal_df = signal_df,
    peak_df = peak_df,
    start_index = start_index,
    end_index = end_index,
    ts_index = ts_index,
    start_pos = start_pos,
    end_pos = end_pos,
    ts_pos = ts_pos,
    direction = read$aln_info$direction,
    filename = filename
  )
  
  return(result)
}

# --------------- Visualization Function 3: Check Reversed Read ----------------
plotReadConversion <- function(read_before, read_after) {
  
  require(ggplot2)
  require(patchwork)
  
  # Функция построения одного графика
  plotOne <- function(read) {
    
    d <- prepareRead(read)
    
    # График
    p <- ggplot(
      d$signal_df,
      aes(x = trace_pos, y = signal, color = nuc)) +
      # chromatogram
      geom_line(linewidth = 0.5) +
      # Alignment rectangle
      annotate("rect", 
               xmin = d$start_pos, 
               xmax = d$end_pos, 
               ymin = -Inf, 
               ymax = Inf, 
               fill = NA, 
               color = "black", 
               linewidth = 0.8) +
      # ts
      geom_vline(xintercept = d$ts_pos, color = "black", linewidth = 0.8) +
      # Позиции пиков
      scale_x_continuous(
        breaks = d$peak_df$peak_pos,
        labels = d$peak_df$nuc,
        expand = expansion(mult = c(0.01, 0.01))
      ) +
      scale_color_manual(values = c(A = "green", C = "blue", G = "orange", T = "red")) +
      labs(
        title = paste0(d$filename, ": ", d$direction),
        x = "Position",
        y = "Signal intensity",
        color = "Nucleotide"
      ) +
      theme_classic() +
      theme(
        axis.text.x = element_text(vjust = 0.5, hjust = 0.5, size = 5),
        plot.title = element_text(face = "bold"), 
        legend.position = "top")
    
    return(p)
  }
  
  # Построение двух графиков
  p_before <- plotOne(read_before)
  p_after <- plotOne(read_after)
  
  # Вывод друг под другом
  (p_before) / (p_after)
}

# ---------------------------- Function 4: Batch Run ---------------------------
signalConversionRun <- function(direction_info, visualization = TRUE) {
  
  control_batch_before <- direction_info$control_info
  sample_batch_before <- direction_info$sample_info
  
  # Reverse conversion
  control_batch_after <- lapply(
    control_batch_before,
    reverseRead
  )
  
  sample_batch_after <- lapply(
    sample_batch_before,
    reverseRead
  )
  
  # Записываем преобразованные reads обратно в direction_info
  direction_info$control_info <- control_batch_after
  direction_info$sample_info <- sample_batch_after
  
  # Визуализация
  if (visualization == TRUE) {
    control_plots <- lapply(
      seq_along(control_batch_before),
      function(i) {
        plotReadConversion(
          read_before = control_batch_before[[i]],
          read_after = control_batch_after[[i]]
        )
      }
    )
    
    sample_plots <- lapply(
      seq_along(sample_batch_before),
      function(i) {
        plotReadConversion(
          read_before = sample_batch_before[[i]],
          read_after = sample_batch_after[[i]]
        )
      }
    )
    
    direction_info$conversion_plots <- list(
      control = control_plots,
      sample = sample_plots
    )
  }
  
  return(direction_info)
}

# ==============================================================================
# -------------------------- Height compensation Step 4 ------------------------
# ==============================================================================

# ----------------------------- Function 1: Mapping ----------------------------
mapping <- function(aln, raw_indices, raw_pos) {
  
  iupac_codes <- list(
    R = c("A", "G"),
    Y = c("C", "T"),
    S = c("G", "C"),
    W = c("A", "T"),
    K = c("G", "T"),
    M = c("A", "C"),
    B = c("C", "G", "T"),
    D = c("A", "G", "T"),
    H = c("A", "C", "T"),
    V = c("A", "C", "G")
  )
  
  called_aligned <- strsplit(as.character(alignedSubject(aln)), "", fixed = TRUE)[[1]]
  raw_aligned <- strsplit(as.character(alignedPattern(aln)), "", fixed = TRUE)[[1]]
  
  called_counter <- 0L
  raw_counter <- 0L
  
  result <- vector("list", length(called_aligned))
  
  for (i in seq_along(called_aligned)) {
    
    nuc_basecalled <- called_aligned[i]
    nuc_raw <- raw_aligned[i]
    
    # Count positions in basecalled sequence
    if (nuc_basecalled != "-") {
      called_counter <- called_counter + 1L
    }
    
    # Count positions in raw sequence
    if (nuc_raw != "-") {
      raw_counter <- raw_counter + 1L
    }
    
    # Nucleotide positions
    nt_basecalled <- if (nuc_basecalled != "-") {
      called_counter
    } else {
      NA_integer_
    }
    
    nt_raw <- if (nuc_raw != "-") {
      raw_indices[raw_counter]
    } else {
      NA_integer_
    }
    
    # Corrected nucleotide
    nuc_corrected <- if (nuc_raw == "-") {
      NA_character_ # gap
    } else if (nuc_basecalled == "-") {
      nuc_raw
    } else if (nuc_raw == "N") {
      nuc_basecalled # N nuc
    } else if (nuc_basecalled %in% iupac_codes[[nuc_raw]]) {
      nuc_basecalled # IUPAC code
    } else if (nuc_raw == nuc_basecalled) {
      nuc_raw # same
    } else if (nuc_raw != nuc_basecalled) {
      nuc_basecalled # mismatch
    } else {
      nuc_raw
    }
    
    # Pseudo-center position
    pseudocenter_pos <- if (!is.na(nt_raw) && nt_raw >= 1 && nt_raw <= length(raw_pos)) {
      raw_pos[nt_raw]
    } else {
      NA_real_
    }
    
    result[[i]] <- data.frame(
      nuc_raw = nuc_raw,
      nuc_basecalled = nuc_basecalled,
      nuc_corrected = nuc_corrected,
      nt_raw = nt_raw,
      pseudocenter_pos = pseudocenter_pos,
      stringsAsFactors = FALSE
    )
  }
  
  result <- do.call(rbind, result)
  rownames(result) <- NULL
  
  return(result)
}

# ----------------------- Function 2: Get pseudo centers -----------------------
getPseudoCenters <- function(read, reference) {
  
  # Input
  raw_seq <- read$main_info$seq
  raw_pos <- read$main_info$pos
  basecalled_seq <- read$main_info$seq_basecalled
  
  aln_raw_ref <- pairwiseAlignment(
    pattern = reference, 
    subject = raw_seq, 
    type = "local")
  
  # Aligned raw sequence
  raw_trimmed <- gsub("-", "", as.character(alignedSubject(aln_raw_ref)))
  raw_start <- start(subject(aln_raw_ref))
  
  raw_indices <- raw_start +
    seq_along(strsplit(raw_trimmed, "", fixed = TRUE)[[1]]) - 1L
  
  # Safety check
  if (any(raw_indices > length(raw_seq))) {
    stop("raw_indices exceed the length of raw_seq.")
  }
  
  aln_called_raw <- pairwiseAlignment(
    pattern = raw_trimmed,
    subject = basecalled_seq,
    type = "local")
  
  mapping <- mapping(
    aln = aln_called_raw,
    raw_indices = raw_indices,
    raw_pos = raw_pos
  )
  
  return(mapping)
}

# ----------------------- Function 3: Get real peak heights -----------------------
getRealPeakHeights <- function(pseudo_centers, sigs_matrix, radius = 4, base_order = c("A", "C", "G", "T")) {
  
  result <- lapply(seq_len(nrow(pseudo_centers)), function(i) {
    
    nuc <- pseudo_centers$nuc_corrected[i]
    #print(paste0("OKAY: ", nuc))
    pos <- pseudo_centers$pseudocenter_pos[i]
    nt_raw <- pseudo_centers$nt_raw[i]
    
    start_cent_idx <- max(1, pos - radius)
    end_cent_idx <- min(nrow(sigs_matrix), pos + radius)
    idx <- start_cent_idx:end_cent_idx
    
    # Получаем сигнал соответствующего канала.
    nuc_idx <- match(nuc, base_order)
    signal <- sigs_matrix[idx, nuc_idx]
    
    # Индекс максимального сигнала в исходной матрице.
    max_idx <- idx[which.max(signal)]
    pseudo_center_sig <- sigs_matrix[pos, nuc_idx]
    real_center_sig <- sigs_matrix[max_idx, nuc_idx]
    
    data.frame(
      nuc_corrected = nuc,
      pseudo_center_pos = pos,
      real_center_pos = max_idx,
      nt_raw = nt_raw,
      diff = max_idx - pos,
      pseudo_center_sig = pseudo_center_sig,
      real_center_sig = real_center_sig
    )
  })
  
  real_heights <- do.call(rbind, result)
  
  return(real_heights)
}


# ------------------------ Function 4: Get best fitting ------------------------
library(minpack.lm)
library(dplyr)

fitPeakDistributionModels <- function(df, 
                                      x_col = "real_center_pos",
                                      y_col = "real_center_sig",
                                      nuc_col = "nuc_corrected", 
                                      max_exp = 4,
                                      n_starts = 30,
                                      delta_bic = 2,
                                      seed = 123
                                      ) {
  
  ################################### - 1 - ####################################
  # Libraries
  require(minpack.lm)
  require(dplyr)
  
  # Constant parameters
  set.seed(seed)
  nucleotides <- c("A", "C", "G", "T")
  
  # Подготовка results
  all_results <- list()
  all_parameters <- list()
  plot_data <- list()
  all_fits <- list()
  selected_models <- list()
  
  ################################### - 2 - ####################################
  
  # Функция для построения модели из k экспонент
  make_formula <- function(k) {
    terms <- paste0("A", 1:k, " * exp(-x / tau", 1:k, ")", collapse = " + ")
    as.formula(paste0("y ~ B + ", terms))
  }
  
  # Функция начальных значений
  make_start <- function(k, y) {
    y_min <- min(y, na.rm = TRUE)
    y_max <- max(y, na.rm = TRUE)
    B_start <- max(0, y_min)
    amplitude <- (y_max - B_start) / k
    starts <- list(B = B_start)
    tau_values <- exp(seq(log(0.02), log(1.5), length.out = k))
    
    for (i in 1:k) {
      starts[[paste0("A", i)]] <- amplitude
      # tau распределяем по разным масштабам
      starts[[paste0("tau", i)]] <- tau_values[i]
    }
    
    return(starts)
  }
  
  ################################### - 3 - ####################################
  
  # --- 1. Обработка каждого нуклеотида ----------------------------------------
  for (nuc in nucleotides) {
    
    ################################# - 3.1 - ##################################
    data_nuc <- df[df[[nuc_col]] == nuc, , drop = FALSE]
    data_nuc <- data_nuc[is.finite(data_nuc[[x_col]]) & is.finite(data_nuc[[y_col]]), , drop = FALSE]
    ################################# - 3.2 - ##################################
    if (nrow(data_nuc) < 10) {
      warning("Not enough observations for nucleotide ", nuc)
      next
    }
    ################################# - 3.3 - ##################################
    # Нормировка x в диапазон 0–1
    x_original <- data_nuc[[x_col]]
    x_min <- min(x_original)
    x_max <- max(x_original)
    data_nuc$x <- (x_original - x_min) / (x_max - x_min)
    data_nuc$y <- data_nuc[[y_col]]
    
    ################################# - 3.4 - ##################################
    # Фиттинг экспоненциальных моделей
    for (k in 1:max_exp) {
      ############################### - 3.4.1 - ################################
      formula_model <- make_formula(k)
      fits <- list()
      y_min <- min(data_nuc$y)
      y_max <- max(data_nuc$y)
      
      ############################### - 3.4.2 - ################################
      # Детерминированный старт
      starts_list <- list(make_start(k, data_nuc$y))
      
      # Несколько случайных стартов
      for (s in 2:n_starts) {
        start <- list(B = runif(1, min = 0, max = max(y_min, y_max * 0.3)))
        
        for (i in seq_len(k)) {
          start[[paste0("A", i)]] <- runif(1, min = -y_max, max = y_max) / k
          start[[paste0("tau", i)]] <- exp(runif(1, log(0.01), log(2)))
          }
        
        starts_list[[length(starts_list) + 1]] <- start
        
        ############################# - 3.4.2.1 - ##############################
        # Bounds
        lower <- list(B = 0)
        upper <- list(B = y_max * 2)
        
        for (i in seq_len(k)) {
          lower[[paste0("A", i)]] <- -y_max * 3
          upper[[paste0("A", i)]] <- y_max * 3
          lower[[paste0("tau", i)]] <- 0.001
          upper[[paste0("tau", i)]] <- 10
        }
        
        ############################# - 3.4.2.2 - ##############################
        # Fit
        fit_try <- try(
          nlsLM(formula_model,
                data = data_nuc,
                start = start,
                lower = unlist(lower),
                upper = unlist(upper),
                control = nls.lm.control(maxiter = 1024)
                ),
          silent = TRUE
          )
        
        ############################# - 3.4.2.2 - ##############################
        if (!inherits(fit_try, "try-error")) {
          fits[[length(fits) + 1]] <- fit_try
        }
      }
      
      ############################### - 3.4.3 - ################################
      # Если ни один fit не получился
      if (length(fits) == 0) {
        warning("No successful fit for ", nuc, ", ", k, " exponent(s)")
        next
      }
      
      ############################### - 3.4.2 - ################################
      # Поиск лучшего fit'а
      rss_values <- sapply(fits, function(fit) {sum(residuals(fit)^2)})
      best_fit <- fits[[which.min(rss_values)]]
      
      x_prediction <- seq(min(data_nuc$x), max(data_nuc$x), length.out = 500)
      prediction_data <- data.frame(
        real_x = x_prediction * (x_max - x_min) + x_min,
        predicted_y = predict(best_fit, newdata = data.frame(x = x_prediction))
      )
      
      ############################### - 3.4.3 - ################################
      # Метрики модели
      residuals_model <- residuals(best_fit)
      
      n <- length(data_nuc$y)
      p <- length(coef(best_fit))
      
      RSS <- sum(residuals_model^2)
      TSS <- sum((data_nuc$y - mean(data_nuc$y))^2)
      R2 <- 1 - RSS / TSS
      RMSE <- sqrt(mean(residuals_model^2))
      Deviation <- sqrt(RSS / (n - p))
      AIC_value <- n * log(RSS / n) + 2 * p
      BIC_value <- n * log(RSS / n) + log(n) * p
      
      ############################### - 3.4.4 - ################################
      plot_data[[paste0(nuc, "_", k)]] <- list(
        nucleotide = nuc,
        model = paste0(k, "_exponential"),
        n_exp = k,
        observed = data.frame(x = data_nuc[[x_col]], y = data_nuc[[y_col]]),
        predicted = prediction_data,
        metrics = data.frame(R2 = R2,
                             RMSE = RMSE,
                             AIC = AIC_value,
                             BIC = BIC_value,
                             Deviation = Deviation)
      )
      
      ############################### - 3.4.5 - ################################
      # Таблица результатов
      all_results[[length(all_results) + 1]] <- data.frame(
        Nucleotide = nuc,
        Model = paste0(k, "_exponential"),
        N = n,
        Parameters = p,
        RSS = RSS,
        RMSE = RMSE,
        R2 = R2,
        AIC = AIC_value,
        BIC = BIC_value,
        Deviation = Deviation
        )
      
      ############################### - 3.4.6 - ################################
      # Параметры модели
      coef_values <- coef(best_fit)
      
      all_parameters[[length(all_parameters) + 1]] <- data.frame(
        Nucleotide = nuc,
        Model = paste0(k, "_exponential"),
        Parameter = names(coef_values),
        Estimate = as.numeric(coef_values)
        )
      
      # Сохраняем fit
      all_fits[[paste0(nuc, "_", k)]] <- list(
        fit = best_fit,
        x_min = x_min,
        x_max = x_max,
        n_exp = k
      )
    }
    
    ################################# - 3.5 - ##################################
    # Выбор модели
    nuc_results <- bind_rows(all_results)
    nuc_results <- nuc_results[nuc_results$Nucleotide == nuc, , drop = FALSE]
    min_BIC <- min(nuc_results$BIC)
    nuc_results$Delta_BIC <- nuc_results$BIC - min_BIC
    
    # Все модели в пределах delta_bic
    acceptable <- nuc_results[nuc_results$Delta_BIC <= delta_bic, ,drop = FALSE]
    
    # Выбираем наиболее простую модель
    acceptable$N_exp <- as.numeric(sub("_exponential", "",acceptable$Model))
    selected_row <- acceptable[which.min(acceptable$N_exp), , drop = FALSE]
    selected_model_name <- selected_row$Model
    selected_n_exp <- selected_row$N_exp
    
    # Сохраняем выбранную модель
    selected_models[[nuc]] <- list(
      nucleotide = nuc,
      model = selected_model_name,
      n_exp = selected_n_exp,
      BIC = selected_row$BIC,
      Delta_BIC = selected_row$Delta_BIC,
      AIC = selected_row$AIC,
      RMSE = selected_row$RMSE,
      R2 = selected_row$R2,
      Deviation = selected_row$Deviation,
      fit = all_fits[[paste0(nuc, "_", selected_n_exp)]]$fit,
      x_min = x_min,
      x_max = x_max
      )
    }
  
  ################################### - 4 - ####################################
  
  # Итоговые таблицы
  results_table <- bind_rows(all_results)
  
  # Добавляем Delta AIC / Delta BIC
  results_table <- results_table %>%
    group_by(Nucleotide) %>%
    mutate(
      Delta_AIC = AIC - min(AIC),
      Delta_BIC = BIC - min(BIC)
    ) %>%
    ungroup()
  
  parameters_table <- bind_rows(all_parameters)
  
  return(
    list(
      results = results_table, # все модели
      parameters = parameters_table, # параметры всех моделей
      fits = all_fits, # все fit-объекты
      plot_data = plot_data, # Данные для визуализации
      selected_models = selected_models # только выбранные модели
    )
  )
}

# ------------------------ Visual Function 5: One model ------------------------
plotPeakModel <- function(plot_data, nucleotide, n_exp, y_max) {
  
  require(ggplot2)
  
  nucleotide_colors <- c(A = "darkgreen", C = "blue", G = "orange", T = "red")
  
  key <- paste0(nucleotide, "_", n_exp)
  data <- plot_data[[key]]
  
  if (is.null(data)) {
    stop("No plot data for ", nucleotide, " with ", n_exp, " exponential(s).")
  }
  
  nuc_color <- nucleotide_colors[nucleotide]
  
  p <- ggplot(
    data$observed,
    aes(x = x, y = y)
    ) +
    geom_point(
      color = nuc_color,
      size = 1.5,
      alpha = 0.5
    ) +
    geom_line(
      data = data$predicted,
      aes(x = real_x, y = predicted_y),
      color = nuc_color,
      linewidth = 1
    ) +
    labs(
      title = paste0(nucleotide, " — ", n_exp, " exponential"),
      subtitle = paste0("BIC = ", round(data$metrics$BIC, 1)),
      x = "Center coordinate",
      y = "Peak height"
    ) +
    
    theme_classic() +
    scale_y_continuous(
      limits = c(0, y_max),
      expand = c(0, 0)
    )
  
  return(p)
}

# ----------------------- Visual Function 6: All models ------------------------
plotAllPeakModels <- function(peak_models, bases = c("A", "C", "G", "T")) {
  
  require(patchwork)
  
  all_y <- c()
  
  for (nuc in bases) {
    for (k in 1:4) {
      key <- paste0(nuc, "_", k)
      plot_info <- peak_models$plot_data[[key]]
      if (!is.null(plot_info)) {
        all_y <- c(all_y, plot_info$observed$y)
      }
    }
  }
  
  y_max <- max(all_y, na.rm = TRUE) * 1.05
  
  plots <- list()
  
  for (nuc in bases) {
    for (k in 1:4) {
      plots[[paste0(nuc, "_", k)]] <- plotPeakModel(peak_models$plot_data,
                                                    nucleotide = nuc,
                                                    n_exp = k,
                                                    y_max = y_max)
    }
  }
  
  combined_plot <- wrap_plots(plots, ncol = 4)
  
  #ggsave("all_peak_models.png", combined_plot, width = 12, height = 10, dpi = 300)
  
  return(combined_plot)
}

# ---------------------- Visual Function 7: Best models ------------------------
plotSelectedPeakModels <- function(peak_models, bases = c("A", "C", "G", "T")) {
  
  require(patchwork)
  
  all_y <- c()
  
  for (nuc in nucleotides) {
    selected <- peak_models$selected_models[[nuc]]
    plot_info <- peak_models$plot_data[[paste0(nuc, "_", selected$n_exp)]]
    all_y <- c(all_y, plot_info$observed$y)
  }
  
  y_max <- max(all_y, na.rm = TRUE) * 1.05
  
  plots <- list()
  
  for (nuc in bases) {
    selected <- peak_models$selected_models[[nuc]]
    n_exp <- selected$n_exp
    plot_info <- peak_models$plot_data[[paste0(nuc, "_", n_exp)]]
    
    plots[[nuc]] <- plotPeakModel(peak_models$plot_data,
                                  nucleotide = nuc,
                                  n_exp = n_exp,
                                  y_max = y_max)
  }
  
  combined_plot <- (plots$A) | (plots$C) | (plots$G) | (plots$T)
  
  #ggsave("selected_peak_models.png", combined_plot, width = 12, height = 3, dpi = 300)
  
  return(combined_plot)

}

################################################################################

library(pracma)
library(tidyr)

# ---------------------- Helper Function 8: Signal factor ----------------------
getRegionIdx <- function(raw_pos, start, end) {
  which(raw_pos >= min(start, end) & raw_pos <= max(start, end))
}

# ---------------------- Helper Function 9: Signal factor ----------------------
getExtendedRegionIdx <- function(raw_pos, pos1, pos2, extension) {
  search_start <- min(pos1, pos2) - extension
  search_end <- max(pos1, pos2) + extension
  
  getRegionIdx(raw_pos = raw_pos, start = search_start, end = search_end)
}

# --------------------- Helper Function 10: Signal factor ----------------------
selectClosestToMidpoint <- function(positions, pos1, pos2) {
  midpoint <- (pos1 + pos2) / 2
  positions[which.min(abs(positions - midpoint))]
}

interpolateCrossing <- function(x1, x2, d1, d2) {
  
  if (d1 == 0) {
    return(x1)
  }
  if (d2 == 0) {
    return(x2)
  }
  
  x1 + (0 - d1) * (x2 - x1) / (d2 - d1)
}

# --------------------- Helper Function 11: Signal factor ----------------------
getCrossingPositions <- function(raw_pos, region_idx, diff, crossing) {
  
  sapply(crossing, function(i) {
    
    x1 <- raw_pos[region_idx[i]]
    x2 <- raw_pos[region_idx[i + 1]]
    
    d1 <- diff[i]
    d2 <- diff[i + 1]
    
    interpolateCrossing(x1 = x1, x2 = x2, d1 = d1, d2 = d2)
  })
}

# --------------------- Helper Function 12: Signal factor ----------------------
findValley <- function(signal, raw_pos, smooth_k = 5) {
  
  if (length(signal) < smooth_k) {
    return(numeric(0))
  }
  
  # Сглаживание
  signal_smooth <- stats::filter(signal, rep(1 / smooth_k, smooth_k), sides = 2)
  signal_smooth <- as.numeric(signal_smooth)
  
  valid <- !is.na(signal_smooth)
  
  if (sum(valid) < 3) {
    return(numeric(0))
  }
  
  valid_idx <- which(valid)
  
  # Не рассматриваем края сглаженного сигнала
  if (length(valid_idx) < 3) {
    return(numeric(0))
  }
  
  # Минимум сглаженного сигнала
  min_idx <- valid_idx[which.min(signal_smooth[valid_idx])]
  
  return(raw_pos[min_idx])
}


# ----- Function 13: Signal factor: find boundary between different Nucs  ------
findBoundaryDiffNucs <- function(sigs_matrix, raw_pos,
                                 nuc1, pos1, nuc2, pos2,
                                 base_order = c("A", "C", "G", "T"),
                                 extension = 5, smooth_k = 3) {
  
  nuc1_idx <- match(nuc1, base_order)
  nuc2_idx <- match(nuc2, base_order)
  
  if (is.na(nuc1_idx) || is.na(nuc2_idx)) {
    return(NA_real_)
  }
  
  # Область строго между центрами пиков
  region_idx <- getRegionIdx(raw_pos, pos1, pos2)
  
  if (length(region_idx) < 2) {
    return(NA_real_)
  }
  
  # 1. Crossing двух сигналов
  sig1 <- sigs_matrix[region_idx, nuc1_idx]
  sig2 <- sigs_matrix[region_idx, nuc2_idx]
  
  diff <- sig1 - sig2
  
  crossing <- which(diff[-length(diff)] * diff[-1] <= 0)
  
  if (length(crossing) > 0) {
    
    crossing_pos <- getCrossingPositions(raw_pos = raw_pos, 
                                         region_idx = region_idx, 
                                         diff = diff, 
                                         crossing = crossing)
    
    return(selectClosestToMidpoint(crossing_pos, pos1, pos2))
  }
  
  # 2. Если crossing нет — ищем valley ТОЛЬКО между центрами
  combined_signal <- sig1 + sig2
  
  valley_pos <- findValley(signal = combined_signal, raw_pos = raw_pos[region_idx], smooth_k = smooth_k)
  
  if (length(valley_pos) > 0) {
    return(valley_pos)
  }
  
  # 3. Если valley не найден — берём точку минимального суммарного сигнала между центрами без выхода за пределы.
  valid <- !is.na(combined_signal)
  
  if (sum(valid) >= 2) {
    
    valid_signal <- combined_signal[valid]
    valid_pos <- raw_pos[region_idx][valid]
    
    return(
      valid_pos[which.min(valid_signal)]
    )
  }
  
  return(NA_real_)
}

# -------- Function 14: Signal factor: find boundary between same Nucs  --------
findBoundarySameNucs <- function(sigs_matrix, raw_pos,
                                 nuc, pos1, pos2, extension = 5,
                                 base_order = c("A", "C", "G", "T"),
                                 smooth_k = 3, tolerance = 0.05) {
  
  nuc_idx <- match(nuc, base_order)
  
  if (is.na(nuc_idx)) {
    return(NA_real_)
  }
  
  # Поиск областей минимума
  findMinima <- function(region_idx) {
    
    if (length(region_idx) < smooth_k) {
      return(numeric(0))
    }
    
    signal <- sigs_matrix[region_idx, nuc_idx]
    
    # Сглаживание
    signal_smooth <- stats::filter(signal, rep(1 / smooth_k, smooth_k), sides = 2)
    signal_smooth <- as.numeric(signal_smooth)
    
    # Убираем NA после smoothing
    valid <- !is.na(signal_smooth)
    
    if (sum(valid) < 3) {
      return(numeric(0))
    }
    
    valid_idx <- which(valid)
    valid_signal <- signal_smooth[valid_idx]
    
    # Минимальное значение
    min_value <- min(valid_signal)
    
    # Допуск вокруг глобального минимума
    threshold <- min_value * (1 + tolerance)
    minima_idx <- valid_idx[signal_smooth[valid_idx] <= threshold]
    
    if (length(minima_idx) == 0) {
      return(numeric(0))
    }
    
    # Объединяем соседние точки в отдельные области
    groups <- split(minima_idx, cumsum(c(TRUE, diff(minima_idx) > 1)))
    
    # Для каждой области выбираем центральную позицию
    minima_pos <- sapply(groups, function(idx) {
      # Центральный индекс области минимума
      center_idx <- idx[ceiling(length(idx) / 2)]
      raw_pos[region_idx[center_idx]]
    })
    
    return(as.numeric(minima_pos))
  }
  
  # 1. Ищем минимум между центрами пиков
  region_idx <- getRegionIdx(raw_pos = raw_pos, start = pos1, end = pos2)
  minima_pos <- findMinima(region_idx)
  
  # Оставляем только точки строго между центрами
  minima_pos <- minima_pos[minima_pos > min(pos1, pos2) & minima_pos < max(pos1, pos2)]
  
  if (length(minima_pos) > 0) {
    # Если найдено несколько областей минимума, выбираем самую близкую к midpoint
    return(selectClosestToMidpoint(positions = minima_pos, pos1 = pos1, pos2 = pos2))
  }
  
  # 2. Минимум между центрами не найден → расширяем область
  region_idx <- getExtendedRegionIdx(raw_pos = raw_pos, pos1 = pos1, pos2 = pos2, extension = extension)
  
  minima_pos <- findMinima(region_idx)
  
  if (length(minima_pos) == 0) {
    return(NA_real_)
  }
  
  # 3. В expanded region выбираем ближайший к midpoint
  return(selectClosestToMidpoint(positions = minima_pos, pos1 = pos1, pos2 = pos2))
}

# -------------- Function 10: Signal factor: run for boundaries ----------------
findAllBoundaries <- function(real_heights, sigs_matrix) {
  
  raw_pos <- seq_len(nrow(sigs_matrix))

  # Расчет границ пиков 
  results <- lapply(seq_len(nrow(real_heights) - 1), function(i) {
    
    nuc1 <- real_heights$nuc_corrected[i]
    pos1 <- real_heights$real_center_pos[i]
    nt_read_1 <- real_heights$nt_raw[i]
    
    nuc2 <- real_heights$nuc_corrected[i + 1]
    pos2 <- real_heights$real_center_pos[i + 1]
    nt_read_2 <- real_heights$nt_raw[i + 1]
    
    # Одинаковые нуклеотиды
    if (nuc1 == nuc2) {
      real_boundary <- findBoundarySameNucs(
        sigs_matrix = sigs_matrix,
        raw_pos = raw_pos,
        nuc = nuc1,
        pos1 = pos1,
        pos2 = pos2
      )
    } else {
      real_boundary <- findBoundaryDiffNucs(
        sigs_matrix = sigs_matrix,
        raw_pos = raw_pos,
        nuc1 = nuc1,
        pos1 = pos1,
        nuc2 = nuc2,
        pos2 = pos2
      )
    }  
    
    data.frame(
      nuc1 = nuc1,
      nt_read_1 = nt_read_1,
      real_pos1 = pos1,
      nuc2 = nuc2,
      nt_read_2 = nt_read_2,
      real_pos2 = pos2,
      real_boundary = real_boundary
    )
  })
  
  boundaries <- do.call(rbind, results)
  
  
  return(boundaries)
}

# --------------- Visualization Function 11*: Check boundaries -----------------
plotChromatogram <- function(heights_window, boundaries, sigs_matrix, filename) {
  
  require(ggplot2)
  
  raw_pos <- nrow(sigs_matrix)
  start_pos <- heights_window[1,"real_center_pos"]
  end_pos <- heights_window[nrow(heights_window), "real_center_pos"]
  raw_pos_window <- start_pos:end_pos
  
  channel_colors <- c(A = "darkgreen", C = "blue", G = "orange", T = "red")
  base_order <- c("A", "C", "G", "T")
  
  peak_labels <- seq_len(nrow(heights_window)) - (nrow(heights_window) %/% 2) - 1
  peak_labels <- ifelse(peak_labels > 0, paste0("+", peak_labels), peak_labels)
  
  x_labels <- data.frame(pos_real = heights_window$real_center_pos, label = peak_labels)
  
  df_long <- rbind(
    data.frame(pos = raw_pos_window, signal = sigs_matrix[raw_pos_window, match("A", base_order)], nuc = "A"),
    data.frame(pos = raw_pos_window, signal = sigs_matrix[raw_pos_window, match("C", base_order)], nuc = "C"),
    data.frame(pos = raw_pos_window, signal = sigs_matrix[raw_pos_window, match("G", base_order)], nuc = "G"),
    data.frame(pos = raw_pos_window, signal = sigs_matrix[raw_pos_window, match("T", base_order)], nuc = "T")
  )
  
  y_max <- max(sigs_matrix[raw_pos_window, ], na.rm = TRUE) * 1.05
  
  p <- ggplot(df_long, aes(x = pos, y = signal, color = nuc)) +
    
    # Реальные центры пиков
    geom_vline(data = heights_window,
               aes(xintercept = real_center_pos),
               linetype = "solid",
               linewidth = 0.4,
               color = "gray50") +
    
    # Границы пиков
    geom_vline(data = boundaries,
               aes(xintercept = real_boundary), 
               linetype = "dashed", 
               linewidth = 0.2, 
               color = "gray30") +
    geom_line(linewidth = 0.7) +
    scale_color_manual(values = channel_colors) +
    
    # Подписи X под центрами пиков 
    scale_x_continuous(breaks = x_labels$pos_real) + #, labels = x_labels$label) +
    scale_x_continuous(breaks = x_labels$pos_real) + # , labels = x_labels$label, expand = c(0, 0)) +
    scale_y_continuous(limits = c(0, y_max), expand = c(0, 0)) +
    theme_classic() +
    theme(
      axis.line = element_line(color = "black", linewidth = 0.25),
      panel.border = element_rect(color = "black", fill = NA, linewidth = 1.75),
      axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
      legend.position = "top"
    ) +
    labs(x = "Window Position", y = "Signal", color = "Nucleotide")
  
  #ggsave(filename = paste0("chromatogram_", filename, ".pdf"), p, device = "pdf", width = 48, height = 3)
   
  return(p)
  
}

# ----------- Function 12:Signal factor: get summary about boundaries ----------
collectPeakInfo <- function(heights_window, boundaries, sigs_matrix, base_order = c("A", "C", "G", "T")) {
  
  heights_trimmed <- heights_window[
    2:(nrow(heights_window) - 1),
    ,
    drop = FALSE
  ]
  
  rounded_boundaries <- round(boundaries$real_boundary)
  
  n_peaks <- nrow(heights_trimmed)
  
  peak_start <- numeric(n_peaks)
  peak_end <- numeric(n_peaks)
  
  peak_start[1] <- rounded_boundaries[1] + 1
  
  for (i in seq_len(n_peaks)) {
    
    boundary_index <- i
    peak_end[i] <- rounded_boundaries[boundary_index + 1]
    
    if (i < n_peaks) {
      peak_start[i + 1] <- peak_end[i] + 1
    }
  }
  
  channel_idx <- match(heights_trimmed$nuc_corrected, base_order)
  result <- lapply(seq_len(n_peaks), function(i) {
    
    # Берём все четыре канала для текущего пика
    peak_data <- sigs_matrix[
      peak_start[i]:peak_end[i],
      ,
      drop = FALSE
    ]
    
    data.frame(
      nuc = heights_trimmed$nuc_corrected[i],
      nt_nuc = heights_trimmed$nt_raw[i],
      real_center_pos = heights_trimmed$real_center_pos[i],
      peak_start = peak_start[i],
      peak_end = peak_end[i],
      A = I(list(peak_data[, 1])),
      C = I(list(peak_data[, 2])),
      G = I(list(peak_data[, 3])),
      T = I(list(peak_data[, 4]))
    )
  })
  
  boundaries_summary <- do.call(rbind, result)
  
  return(boundaries_summary)
}

# ----------- Function 13: Signal factor: predict peak distribution ------------
predictPeakDistribution <- function(model_set, nucleotide, x) {
  
  model_info <- model_set[[nucleotide]]
  
  if (is.null(model_info)) {
    stop("No model available for nucleotide ", nucleotide)
  }
  
  # Переводим реальную позицию в 0–1
  x_norm <- (x - model_info$x_min) / (model_info$x_max - model_info$x_min)
  
  # Создаём data.frame для predict()
  newdata <- data.frame(x = x_norm)
  
  # Предсказанная интенсивность
  y_pred <- predict(model_info$fit, newdata = newdata)
  
  return(as.numeric(y_pred))
}

# --------------- Function 14: Signal factor: correction factor ----------------
calculatePosCorrectionFactor <- function(model_set, nucleotide, x, reference_position = NULL) {
  
  if (is.null(reference_position)) {
    reference_position <- model_set[[nucleotide]]$x_min
  }
  
  predicted_signal <- predictPeakDistribution(model_set = model_set,
                                              nucleotide = nucleotide,
                                              x = x)
  
  reference_signal <- predictPeakDistribution(model_set = model_set,
                                              nucleotide = nucleotide,
                                              x = reference_position)
  
  correction_factor <- reference_signal / predicted_signal
  
  return(as.numeric(correction_factor))
}

# --------------- Function 15: Signal factor: signal correction ----------------
signalPosCorrection <- function(peak_info, model_set, nucleotide) {
  
  result <- peak_info
  
  result$signal_corrected <- lapply(
    seq_len(nrow(peak_info)),
    function(i) {
      
      # Реальные позиции всех точек текущего пика
      x <- peak_info$peak_start[i]:peak_info$peak_end[i]
      
      # Исходный сигнал
      signal <- peak_info[[nucleotide]][[i]]
      
      # Correction factor для каждой точки
      correction_factor <- calculatePosCorrectionFactor(
        model_set = model_set,
        nucleotide = nucleotide,
        x = x
      )
      
      # Корректируем весь сигнал
      signal_corrected <- signal * correction_factor
      
      return(signal_corrected)
    }
  )
  
  return(result)
}

# ------------- Function 16: Signal factor: signal correction run --------------
signalPosCorrectionRun <- function(peak_info, model_set) {
  
  base_order <- c("A", "C", "G", "T")
  
  result <- lapply(base_order, function(nucleotide) {
    signalPosCorrection(
        peak_info = peak_info[c("nuc", "nt_nuc", "real_center_pos", "peak_start", "peak_end", nucleotide)],
        model_set = model_set,
        nucleotide = nucleotide
      )
    }
  )
  
  names(result) <- base_order
  
  # Общая информация о пиках
  output <- result[[1]][
    c("nuc", "nt_nuc", "real_center_pos", "peak_start", "peak_end")
  ]
  
  # Добавляем raw и corrected signal для каждого нуклеотида
  for (nucleotide in base_order) {
    
    current <- result[[nucleotide]]
    
    output[[paste0(nucleotide, "_raw_sigs")]] <- 
      current[[nucleotide]]
    
    output[[paste0(nucleotide, "_pos_cor_sigs")]] <- 
      current$signal_corrected
  }
  
  return(output)
}

# ---------- Function 17: Signal factor: signal dye correction factor ----------
getCorrectedPeakHeights <- function(signals_corrected) {
  
  result <- data.frame(
    nuc = signals_corrected$nuc,
    A_height_cor = numeric(nrow(signals_corrected)),
    C_height_cor = numeric(nrow(signals_corrected)),
    G_height_cor = numeric(nrow(signals_corrected)),
    T_height_cor = numeric(nrow(signals_corrected))
  )
  
  bases <- c("A", "C", "G", "T")
  
  for (i in seq_len(nrow(signals_corrected))) {
    
    start <- signals_corrected$peak_start[i]
    center <- signals_corrected$real_center_pos[i]
    
    center_index <- center - start + 1
    
    for (base in bases) {
      
      signal <- signals_corrected[[paste0(base, "_pos_cor_sigs")]][[i]]
      result[i, paste0(base, "_height_cor")] <- signal[center_index]
    }
  }
  
  return(result)
}

# ---------- Function 18: Signal factor: signal dye correction factor ----------
signalDyeCorrectionFactor <- function(pos_corrected_peak_heights) {
  
  bases <- c("A", "C", "G", "T")
  
  # Типичная высота основного пика для каждого нуклеотида
  median_heights <- sapply(bases, function(base) {
    
    values <- pos_corrected_peak_heights[
      pos_corrected_peak_heights$nuc == base,
      paste0(base, "_height_cor")
    ]
    
    values <- values[is.finite(values) & values > 0]
    
    if (length(values) == 0) {
      return(NA_real_)
    }
    
    median(values)
  })
  
  # Общий уровень интенсивности
  reference_height <- exp(mean(log(median_heights), na.rm = TRUE))
  
  # Correction factor
  correction_factors <- reference_height / median_heights
  
  result <- data.frame(
    nuc = bases,
    median_height = as.numeric(median_heights),
    dye_correction_factor = as.numeric(correction_factors)
  )
  
  return(result)
}

# ---------- Function 19: Signal factor: apply_dye_correction ----------
applyDyeCorrectionFactor <- function(signals_pos_corrected,
                                     dye_correction_factors) {
  
  bases <- c("A", "C", "G", "T")
  
  signals_pos_dye_corrected <- signals_pos_corrected
  
  for (base in bases) {
    
    # Получаем correction factor для текущего канала
    correction_factor <- dye_correction_factors[
      dye_correction_factors$nuc == base,
      "dye_correction_factor"
    ]
    
    if (length(correction_factor) != 1 || is.na(correction_factor)) {
      stop(
        paste0(
          "Dye correction factor для канала ",
          base,
          " отсутствует или некорректен."
        )
      )
    }
    
    # Названия колонок
    input_col <- paste0(base, "_pos_cor_sigs")
    output_col <- paste0(base, "_pos_dye_cor_sigs")
    
    # Применяем correction factor ко всем peak windows
    signals_pos_dye_corrected[[output_col]] <- lapply(
      signals_pos_dye_corrected[[input_col]],
      function(signal) {
        signal * correction_factor
      }
    )
  }
  
  return(signals_pos_dye_corrected)
}

# ------------ Function 20: Signal factor: signal dye correction run -----------

signalDyeCorrectionRun <- function(pos_corrected_peak_heights, signals_pos_corrected) {
  
  dye_correction_factors <- signalDyeCorrectionFactor(
    pos_corrected_peak_heights = pos_corrected_peak_heights
    )
  
  signals_pos_dye_corrected <- applyDyeCorrectionFactor(
    signals_pos_corrected = signals_pos_corrected,
    dye_correction_factors = dye_correction_factors
  )
  
  return(signals_pos_dye_corrected)
}


plotRealHeights <- function(real_heights, filename) {
  
  nuc_colors <- c(
    A = "green4",
    C = "blue3",
    G = "orange2",
    T = "red3"
  )
  
  plots <- lapply(c("A", "C", "G", "T"), function(nuc) {
    
    plot_data <- real_heights[
      real_heights$nuc_corrected == nuc,
    ]
    
    ggplot(
      plot_data,
      aes(
        x = nt_raw,
        y = real_center_sig
      )
    ) +
      geom_line(
        linetype = "dashed",
        color = nuc_colors[[nuc]],
        linewidth = 0.45
      ) +
      geom_point(
        color = nuc_colors[[nuc]],
        size = 1.8
      ) +
      labs(
        title = nuc,
        x = "Center position",
        y = "Peak height"
      ) +
      scale_x_continuous(
        expand = expansion(mult = c(0.01, 0.01))
      ) +
      scale_y_continuous(
        limits = c(0, 2500)
      ) +
      theme_bw(base_size = 10) +
      theme(
        plot.title = element_text(
          face = "bold",
          color = nuc_colors[[nuc]]
        )
      )
  })
  
  comdined_plot <- wrap_plots(plots, ncol = 1)
  
  ggsave(filename = paste0(filename, "_height_distribution.png"),
    plot = comdined_plot,
    width = 12,
    height = 8,
    dpi = 300
  )
  
}


###############################################################################
signalPosDyeCorrectionRun <- function(batch,      # set of reads
                                      reference   # DNAString
                                      ) {
  
  results <- list()
  
  for (read in batch) {
    
    filename <- read$main_info$filename
    
    # ---- I ----
    pseudo_centers <- getPseudoCenters(read = read, 
                                       reference = reference)
    #print(tail(pseudo_centers, 200))
    
    # ---- II ----
    real_peak_heights <- getRealPeakHeights(pseudo_centers = pseudo_centers,
                                            sigs_matrix = read$main_info$sigs)
    
    #plot_heights <- plotRealHeights(real_peak_heights, filename)
    
    # ---- III ----
    boundaries <- findAllBoundaries(real_heights = real_peak_heights, 
                                    sigs_matrix = read$main_info$sigs)

    plot_boudaries <- plotChromatogram(heights_window = real_peak_heights,
                                       boundaries = boundaries,
                                       sigs_matrix = read$main_info$sigs,
                                       filename = read$main_info$filename)
    
    #print(tail(boundaries, 20))
    #print(boundaries[!complete.cases(boundaries), ])
    #print(read$main_info$sigs[7141:7170,1:4])
    
    # ---- IV ----
    peak_info <- collectPeakInfo(heights_window = real_peak_heights, 
                                 boundaries = boundaries, 
                                 sigs_matrix = read$main_info$sigs)
    
    
    # ---- [V] ----
    peak_models <- fitPeakDistributionModels(df = real_peak_heights)
    
    #plot_all_models <- plotAllPeakModels(peak_models)
    #plot_selected_models <- plotSelectedPeakModels(peak_models)
    
    #plot_heights_distribution <- signal_correction_factor()
    
    # ---- [VI] ----
    signals_pos_corrected <- signalPosCorrectionRun(peak_info = peak_info, 
                                                    model_set = peak_models$selected_models)
    
    #plot_signal_pos_corrected <- plot_signal_correction_run(signal_pos_corrected)
    
    # ---- {VII} ----
    pos_corrected_peak_heights <- getCorrectedPeakHeights(signals_corrected = signals_pos_corrected)
    
    # ---- [[VIII]] ----
    signals_pos_dye_corrected <- signalDyeCorrectionRun(pos_corrected_peak_heights = pos_corrected_peak_heights, 
                                                        signals_pos_corrected = signals_pos_corrected)
    
    
    #plot_signal_pos_dye_corrected <- Нужна функция!!!
    
    plots <- list(
    #  plot_all_models = plot_all_models,
    #  plot_selected_models = plot_selected_models,
    #  plot_boudaries = plot_boudaries,
    #  plot_heights_distribution = plot_heights_distribution,
    #  plot_signal_pos_corrected = plot_signal_pos_corrected
      #plot_signal_pos_dye_corrected = plot_signal_pos_dye_corrected
    )
    
    results[[filename]] <- list(
      signals_pos_dye_corrected = signals_pos_dye_corrected,
      plots = plots
    )
  }
  
  return(results)
}

################################################################################
getRunInfo <- function(read) {
  
  data.frame(
    direction = as.character(read$aln_info$direction),
    user = as.character(read$additional_info$user),
    machine = as.character(read$additional_info$machine),
    dye_kit = as.character(read$additional_info$dye_kit),
    run_start = as.character(read$additional_info$run_info$start),
    run_stop = as.character(read$additional_info$run_info$stop),
    run_duration = as.character(read$additional_info$run_info$duration),
    stringsAsFactors = FALSE
  )
}

getBatchRunInfo <- function(batch) {
  
  result <- lapply(batch, getRunInfo)
  result <- do.call(rbind, result)
  rownames(result) <- names(batch)
  
  return(result)
}

runCheck <- function(imported) {
  
  # Объединяем control и sample в один список
  batch <- c(
    imported$control_info,
    imported$sample_info
  )
  
  # Получаем информацию о запуске для каждого прочтения
  run_info <- getBatchRunInfo(batch)
  
  # Добавляем название прочтения
  run_info$read <- rownames(run_info)
  rownames(run_info) <- NULL
  
  # Указываем тип прочтения
  run_info$type <- c(
    rep("control", length(imported$control_info)),
    rep("sample", length(imported$sample_info))
  )
  
  # Определяем уникальный sequencing run
  run_id <- interaction(
    run_info$user,
    run_info$machine,
    run_info$dye_kit,
    run_info$run_start,
    run_info$run_stop,
    run_info$run_duration,
    drop = TRUE
  )
  
  # Присваиваем каждому уникальному запуску run 1, run 2, ...
  run_info$run <- paste0(
    "run ",
    match(run_id, unique(run_id))
  )
  
  return(run_info)
}


# 2. Получаем таблицы pos_dye_corrected
getReadData <- function(type, read, pos_dye_corrected, imported) {
  
  if (type == "control") {
    pos_data <- pos_dye_corrected$control[[read]]$signals_pos_dye_corrected
    imported_read <- imported$control_info[[read]]
  } else if (type == "sample") {
    pos_data <- pos_dye_corrected$sample[[read]]$signals_pos_dye_corrected
    imported_read <- imported$sample_info[[read]]
  } else {
    stop("Неизвестный тип прочтения: ", type)
  }
  
  if (is.null(pos_data)) {
    stop("Не найдены данные pos_dye_corrected для: ", read)
  }
  
  if (is.null(imported_read)) {
    stop("Не найдены данные imported для: ", read)
  }
  
  # Проверяем наличие nuc
  if (!"nuc" %in% colnames(pos_data)) {
    stop("В данных ", read, " отсутствует столбец 'nuc'.")
  }
  
  sequence <-  DNAString(paste0(pos_data$nuc, collapse = ""))
  
  # Добавляем nt_nuc.
  # Если он уже есть, используем существующий.
  if (!"nt_nuc" %in% colnames(pos_data)) {
    pos_data$nt_nuc <- seq_len(nrow(pos_data))
  }
  
  # Качество берём из main_info$qual по nt_nuc
  qual <- imported_read$main_info$qual
  
  pos_data$qual <- qual[pos_data$nt_nuc]
  
  # Сохраняем target site
  target_site <- imported_read$aln_info$ts
  
  list(
    read = read,
    type = type,
    sequence = sequence,
    qual = qual,
    target_site = target_site,
    pos_data = pos_data,
    imported_read = imported_read
  )
}

# 6. Выравниваем каждое прочтение относительно reference
makeAlignment <- function(reference, read) {
  
  alignment <- pairwiseAlignment(
    pattern = reference$sequence,
    subject = read$sequence,
    type = "local"
  )
  
  #writePairwiseAlignments(alignment)
  
  aligned_reference <- strsplit(as.character(alignedPattern(alignment)), "")[[1]]
  aligned_read <-strsplit( as.character(alignedSubject(alignment)), "")[[1]]
  
  # Positions
  ref_pos <- integer(length(aligned_reference))
  read_pos <- integer(length(aligned_read))
  
  ref_counter <- 0
  read_counter <- 0
  
  for (i in seq_along(aligned_reference)) {
    
    if (aligned_reference[i] != "-") {
      ref_counter <- ref_counter + 1
      ref_pos[i] <- ref_counter
    }
    
    if (aligned_read[i] != "-") {
      read_counter <- read_counter + 1
      read_pos[i] <- read_counter
    }
  }
  
  # Оставляем только позиции, присутствующие одновременно в reference и прочтении
  keep <- aligned_reference != "-" & aligned_read != "-"
  
  ref_shift <- reference$pos_data$nt_nuc[1] + start(pattern(alignment))
  read_shift <- read$pos_data$nt_nuc[1] + start(subject(alignment))
  
  alignment_data <- data.frame(
    nt_ref = ref_pos[keep] + ref_shift - 2,
    nt_nuc = read_pos[keep] + read_shift - 2,
    nuc = aligned_read[keep],
    stringsAsFactors = FALSE
  )
  
  # Добавляем качество
  alignment_data$qual <- read$qual[alignment_data$nt_nuc]

  alignment_data
}

getPeakHeights <- function(pos_data) {
  
  result <- data.frame(
    nuc = pos_data$nuc,
    nt_nuc = pos_data$nt_nuc,
    A_height = NA_real_,
    C_height = NA_real_,
    G_height = NA_real_,
    T_height = NA_real_
  )
  
  bases <- c("A", "C", "G", "T")
  
  for (i in seq_len(nrow(pos_data))) {
    
    peak_positions <- pos_data$peak_start[i]:pos_data$peak_end[i]
    center_index <- which(peak_positions == pos_data$real_center_pos[i])
    
    if (length(center_index) != 1) {
      next
    }
    
    for (base in bases) {
      signal <- pos_data[[paste0(base, "_pos_dye_cor_sigs")]][i]
      result[[paste0(base, "_height")]][i] <- signal[[1]][[center_index]]
    }
  }
  
  return(result)
}

getAlignedPeakHeights <- function(alignment_data, pos_data, mode) {
  
  heights <- getPeakHeights(pos_data)
  
  result <- alignment_data[, mode, drop = FALSE]
  colnames(result) <- "nt_nuc"
  result <- merge(result, heights, by = "nt_nuc", all.x = TRUE, sort = FALSE)
  
  return(result)
}

prepareRunCorrectionData <- function(run_info, imported, pos_dye_corrected, 
                                     quality_threshold = 20) {
  
  # 1. Выбираем первое прочтение из run 1 в качестве референса
  run1 <- run_info[run_info$run == "run 1", , drop = FALSE]
  reference_info <- run1[1, , drop = FALSE]
  
  reference <- getReadData(
    type = reference_info$type, 
    read = reference_info$read, 
    pos_dye_corrected = pos_dye_corrected, 
    imported = imported)
  
  # 2. Выбираем остальные прочтения с отличным run
  other_info <- run_info[run_info$run != "run 1", , drop = FALSE]
  
  other_reads <- lapply(seq_len(nrow(other_info)), function(i) {
    getReadData(
      type = other_info$type[i], 
      read = other_info$read[i], 
      pos_dye_corrected = pos_dye_corrected, 
      imported = imported)
  })
  
  # 3. Выравнивания
  alignmets <- lapply(other_reads, function(read) {
    makeAlignment(
      reference = reference,
      read = read
    )
  })
  
  names(alignmets) <- other_info$read
  
  # 4. Объединение выравниваний
  alignments_named <- lapply(seq_along(alignmets), function(i) {
    data <- alignmets[[i]]
    data <- data[, c("nt_ref", "nt_nuc", "nuc", "qual"), drop = FALSE]
    colnames(data) <- c("nt_ref",
                        paste0("nt_nuc_read", i),
                        paste0("nuc_read", i),
                        paste0("qual_read", i)
                        )
    return(data)
  })
  
  merged_alignments <- Reduce(
    function(x, y) {
      merge(x, y, by = "nt_ref", all = FALSE, sort = TRUE)
    },
    alignments_named
  )
  
  rownames(merged_alignments) <- NULL
  
  nuc_cols <- grep("^nuc_", colnames(merged_alignments),  value = TRUE)
  nt_nuc_cols <- grep("^nt_nuc_", colnames(merged_alignments),  value = TRUE)
  qual_cols <- grep("^qual_", colnames(merged_alignments),  value = TRUE)
  reference_nuc <- nuc_cols[1]
  
  merged_alignments <- merged_alignments[, c(
    reference_nuc,
    "nt_ref",
    nt_nuc_cols,
    qual_cols
  ), drop = FALSE]
  
  merged_alignments$qual_ref <- reference$qual[merged_alignments$nt_ref]
  merged_alignments <- merged_alignments %>% rename(nuc_aln = nuc_read1)
  
  # 5. Фильтрация
  merged_alignments_filtered <- merged_alignments[merged_alignments$nt_ref != reference$target_site,]
  
  qual_cols <- grep("^qual_read", colnames(merged_alignments_filtered), value = TRUE)
  merged_alignments_filtered <- merged_alignments_filtered[
    merged_alignments_filtered$qual_ref >= quality_threshold &
      apply(
        merged_alignments_filtered[, qual_cols, drop = FALSE],
        1,
        function(x) all(x >= quality_threshold)
      ),
  ]
  
  # 6. Get Peak Heights
  height_tables <- list()
  
  height_tables$ref <- getAlignedPeakHeights(
    alignment_data = merged_alignments_filtered,
    pos_data = reference$pos_data,
    mode = "nt_ref"
  )
  
  for (i in seq_along(other_reads)) {
    read_name <- paste0("read", i)
    height_tables[[read_name]] <- getAlignedPeakHeights(
      alignment_data = merged_alignments_filtered,
      pos_data = other_reads[[i]]$pos_data,
      mode = paste0("nt_nuc_", read_name)
    )
  }
  
  return(height_tables)
}

getMajorHeights <- function(ref_table, read_table) {
  
  bases <- c("A", "C", "G", "T")
  height_cols <- paste0(bases, "_height")
  
  major_base <- apply(
    ref_table[, height_cols, drop = FALSE],
    1,
    function(x) {
      bases[which.max(x)]
    }
  )
  
  read_major_height <- mapply(
    function(base, i) {
      read_table[[paste0(base, "_height")]][i]
    },
    major_base,
    seq_len(nrow(read_table))
  )
  
  result <- data.frame(
    nt_aln = ref_table$nt_nuc,
    major_base = major_base,
    ref_height = pmax(
      ref_table$A_height,
      ref_table$C_height,
      ref_table$G_height,
      ref_table$T_height,
      na.rm = TRUE
    ),
    read_height = read_major_height
  )
  
  return(result)
}

signalRunCorrectionFactor <- function(prepared_data) {
  
  srcf <- list()
  
  for (i in 2:length(prepared_data)) {
    
    read_name <- names(prepared_data)[i]
    
    major_data <- getMajorHeights(
      ref_table = prepared_data$ref,
      read_table = prepared_data[[read_name]]
    )
    
    major_data$ratio <- (
      major_data$ref_height /
        major_data$read_height
    )
    
    srcf[[read_name]] <- median(
      major_data$ratio,
      na.rm = TRUE
    )
  }
  
  return(srcf)
  
}

applySRCF <- function(run_info, pos_dye_corrected, run_correction_factors) {
  
  bases <- c("A", "C", "G", "T")
  
  other_info <- run_info[run_info$run != "run 1", , drop = FALSE]
  
  result <- list()
  
  for (i in (1:nrow(other_info))) {
    read_name <- other_info$read[i]
    read_type <- other_info$type[i]
    pos_data <-pos_dye_corrected[[read_type]][[read_name]]$signals_pos_dye_corrected
    srcf <- run_correction_factors[[i]]
    
    for (base in bases) {
      signal_col <- paste0(base, "_pos_dye_cor_sigs")
      corrected_col <- paste0(base, "_pos_dye_run_cor_sigs")
      pos_data[[corrected_col]] <- lapply(
        pos_data[[signal_col]], 
        function(signal) {
          signal * srcf
        }
      )
    }
    pos_dye_corrected[[read_type]][[read_name]]$signals_pos_dye_corrected <- pos_data
  }
  
  return(pos_dye_corrected)
}

dataFormation <- function(pos_dye_run_corrected) {
  
  bases <- c("A", "C", "G", "T")
  common_cols <- c("nuc", "nt_nuc", "real_center_pos", "peak_start", "peak_end")
  corrected <- list()
  
  for (type in c("control", "sample")) {
    
    corrected[[type]] <- list()
    reads <- names(pos_dye_run_corrected[[type]])
    
    for (read in reads) {
      
      # Исходная таблица для данного read
      signals <- pos_dye_run_corrected[[type]][[read]]$signals_pos_dye_corrected
      
      # 1. RAW SIGNALS
      raw_cols <- paste0(bases, "_raw_sigs")
      raw_table <- signals[, c(common_cols, raw_cols), drop = FALSE]
      colnames(raw_table) <- c(common_cols, paste0(bases, "_sigs"))
      
      # 2. POSITION CORRECTED SIGNALS
      pos_cols <- paste0(bases, "_pos_cor_sigs")
      pos_table <- signals[, c(common_cols, pos_cols), drop = FALSE]
      colnames(pos_table) <- c(common_cols, paste0(bases, "_sigs"))
      
      # 3. POSITION + DYE CORRECTED SIGNALS
      pos_dye_cols <- paste0(bases, "_pos_dye_cor_sigs")
      pos_dye_table <- signals[, c(common_cols, pos_dye_cols), drop = FALSE]
      colnames(pos_dye_table) <- c(common_cols, paste0(bases, "_sigs"))
      
      # 4. POSITION + DYE + RUN CORRECTED SIGNALS
      run_cols <- paste0(bases, "_pos_dye_run_cor_sigs")
      
      # Проверяем наличие всех run-corrected столбцов
      has_run_correction <- all(run_cols %in% colnames(signals))
      
      if (has_run_correction) {
        
        run_table <- signals[, c(common_cols, run_cols), drop = FALSE]
        colnames(run_table) <- c(common_cols, paste0(bases, "_sigs"))
        status <- "run"
        
      } else {
        run_table <- NULL
        status <- "pos_dye"
      }
      
      # Сохраняем результат для данного read
      corrected[[type]][[read]] <- list(
        raw_signals = raw_table,
        pos_corrected_signals = pos_table,
        pos_dye_corrected_signals = pos_dye_table,
        pos_dye_run_corrected_signals = run_table,
        status = status
      )
    }
  }
  
  return(corrected)
}

signalRunCorrectionRun <- function(conversed_info, pos_dye_corrected, run_info) {
  
  prepared_data <- prepareRunCorrectionData(
    run_info = run_info, 
    imported = conversed_info, 
    pos_dye_corrected = pos_dye_corrected,
    quality_threshold = 20
    )
  
  run_correction_factors <- signalRunCorrectionFactor(prepared_data)
  
  pos_dye_run_corrected <- applySRCF(
    run_info = run_info, 
    pos_dye_corrected = pos_dye_corrected, 
    run_correction_factors = run_correction_factors)
  
  corrected <- dataFormation(pos_dye_run_corrected)
  
  return(corrected)
}





# ------------------- Function 23: Signal Correction Main Run ------------------
signalCorrectionRun <- function(conversed_info, visualization = TRUE) {
  
  control_batch <- conversed_info$control_info
  sample_batch <- conversed_info$sample_info
  reference <- conversed_info$reference_target_info$reference_info
  
  # --- 1. Position & Dye Correction -------------------------------------------
  control_pos_dye_corrected_signals <- signalPosDyeCorrectionRun(
    batch = control_batch,
    reference = reference)
  
  sample_pos_dye_corrected_signals <- signalPosDyeCorrectionRun(
    batch = sample_batch,
    reference = reference)
  
  pos_dye_corrected <- list(
    control = control_pos_dye_corrected_signals,
    sample = sample_pos_dye_corrected_signals,
    status = "pos_dye_corrected"
  )
  
  # --- 2. Run Correction (if needed) ------------------------------------------
  run_info <- runCheck(imported = conversed_info)

  if (length(unique(run_info$run)) > 1) {
    
    pos_dye_run_corrected <- signalRunCorrectionRun(
      conversed_info = conversed_info,
      pos_dye_corrected = pos_dye_corrected,
      run_info = run_info
      )
    
    #plotCorrectedSignals(corrected = pos_dye_run_corrected,
    #                     type = "sample",
    #                     read_name = names(pos_dye_run_corrected$sample)[[2]],
    #                     filename = "sample_pos_dye_run_plot.png")
    
    
    #plotCorrectedSignals(corrected = pos_dye_run_corrected,
    #                     type = "control",
    #                     read_name = names(pos_dye_run_corrected$control)[[1]],
    #                     filename = "control_pos_dye_run_plot.png")
    
    return(pos_dye_run_corrected)
    
  } else {
    
    pos_dye_corrected <- dataFormation(pos_dye_corrected)
    
    #plotCorrectedSignals(corrected = pos_dye_corrected,
    #                     filename = "pos_dye_plot.png")
    
    return(pos_dye_corrected)
    
  }
}

###############################################################################

# -------- Visualization Function 17: Signal factor: signal correction ---------
plotCorrectedSignals <- function(
    corrected,
    type = c("control", "sample"),
    read_name = NULL,
    filename = NULL,
    width = 8,
    height = 8,
    dpi = 300
) {
  
  library(ggplot2)
  library(patchwork)
  
  type <- match.arg(type)
  
  if (is.null(read_name)) {
    read_name <- names(corrected[[type]])[1]
  }
  
  if (!read_name %in% names(corrected[[type]])) {
    stop("Прочтение '", read_name, "' не найдено в corrected$", type)
  }
  
  read_data <- corrected[[type]][[read_name]]
  
  nuc_colors <- c(A = "green4", C = "blue3", G = "orange2", T = "red3")
  
  line_types <- c(raw = "solid", pos = "dashed", pos_dye = "dotdash", pos_dye_run = "dotted")
  
  collectSignals <- function(read_data, nuc, n = NULL) {
    
    stages <- list(
      raw = read_data$raw_signals,
      pos = read_data$pos_corrected_signals,
      pos_dye = read_data$pos_dye_corrected_signals,
      pos_dye_run = read_data$pos_dye_run_corrected_signals
    )
    
    result <- list()
    
    for (stage in names(stages)) {
      
      tab <- stages[[stage]]
      signal_col <- paste0(nuc, "_sigs")
      
      if (is.null(tab) || !signal_col %in% names(tab)) {
        next
      }
      
      # Извлекаем первые n строк таблицы
      if (!is.null(n)) {
        tab <- tab[488:511, , drop = FALSE]
      }
      
      # Собираем вектор сигналов из выбранных строк
      signal <- unlist(tab[[signal_col]], use.names = FALSE)
      signal <- as.numeric(signal)
      
      if (length(signal) == 0) {
        next
      }
      
      result[[stage]] <- data.frame(
        position = seq_along(signal),
        signal = signal,
        stage = stage
      )
    }
    
    if (length(result) == 0) {
      return(NULL)
    }
    
    do.call(rbind, result)
  }
  
  plots <- list()
  
  for (nuc in c("A", "C", "G", "T")) {
    
    plot_data <- collectSignals(read_data = read_data, nuc = nuc, n = 50)
    
    if (is.null(plot_data)) {
      plots[[nuc]] <- ggplot() +
        labs(title = nuc, subtitle = "Нет данных") +
        theme_void()
      next
    }
    
    # Цвет заливки: серый для raw, цвет нуклеотида для коррекций
    plot_data$fill_group <- ifelse(plot_data$stage == "raw", "raw", nuc)
    
    fill_colors <- c(raw = "grey60", setNames(nuc_colors[[nuc]], nuc))
    
    plots[[nuc]] <- ggplot(
      plot_data,
      aes(x = position, y = signal, group = stage, color = stage, linetype = stage)
      ) +
      geom_ribbon(
        aes(ymin = 0, ymax = signal, fill = fill_group, group = stage),
        alpha = 0.12,
        colour = NA
        ) +
      geom_line(linewidth = 0.45, na.rm = TRUE) +
      scale_color_manual(
        name = "Correction stage",
        values = c(raw = "grey60", pos = nuc_colors[[nuc]], pos_dye = nuc_colors[[nuc]], pos_dye_run = nuc_colors[[nuc]]),
        breaks = names(line_types),
        labels = c(raw = "Raw signal", pos = "SPCF", pos_dye = "SPCF + SDCF", pos_dye_run = "SPCF + SDCF + SRCF"),
        drop = FALSE,
        na.translate = FALSE
      ) +
      scale_linetype_manual(
        name = "Correction stage",
        values = line_types,
        breaks = names(line_types),
        labels = c(raw = "Raw signal", pos = "SPCF", pos_dye = "SPCF + SDCF", pos_dye_run = "SPCF + SDCF + SRCF"),
        drop = FALSE,
        na.translate = FALSE
      ) +
      scale_fill_manual(values = fill_colors, guide = "none") +
      scale_x_continuous(expand = c(0, 0)) +
      scale_y_continuous(limits = c(0, 2000), expand = c(0, 0)) +
      coord_cartesian(expand = FALSE) +
      labs(title = paste0(type," ",nuc), x = "Signal coordinate", y = "Signal intensity") +
      theme_bw(base_size = 10) +
      theme(plot.title = element_text(face = "bold", colour = nuc_colors[[nuc]]))
  }
  
  final_plot <- wrap_plots(plots, ncol = 1)

  if (!is.null(filename)) {
    ggsave(filename = filename, plot = final_plot, width = width, height = height, dpi = dpi)
  }
  
  return(final_plot)
}



# ==============================================================================
# -------------------------- Background Threshold Step -------------------------
# ==============================================================================

# ------------- Function 1: Background threshold:  -------------
getTargetWindow <- function(signals_pos_dye_corrected, read, window = 50) {
  
  # Нужные столбцы
  columns <- c(
    "nuc",
    "nt_nuc",
    "real_center_pos",
    "peak_start",
    "peak_end",
    "A_sigs",
    "C_sigs",
    "G_sigs",
    "T_sigs"
  )
  
  # Проверяем наличие столбцов
  missing_columns <- setdiff(columns, names(signals_pos_dye_corrected))
  
  if (length(missing_columns) > 0) {
    stop(paste("Отсутствуют столбцы:", paste(missing_columns, collapse = ", ")))
  }
  
  # Оставляем только нужные столбцы
  result <- signals_pos_dye_corrected[, columns, drop = FALSE]
  
  # Target site
  target_site <- read$aln_info$ts
  
  if (is.null(target_site) || is.na(target_site)) {
    stop("read$aln_info$ts не задан.")
  }
  
  # Ищем строку, соответствующую Target_site
  target_idx <- which(result$nt_nuc == target_site)
  
  if (length(target_idx) == 0) {
    stop(paste0("Target_site = ", target_site, " не найден в nt_raw."))
  }
  
  # Если почему-то найдено несколько строк
  if (length(target_idx) > 1) {
    warning(paste0("Target_site = ", target_site, " найден в нескольких строках. ", "Используется первая."))
    
    target_idx <- target_idx[1]
  }
  
  # Индексы window строк до и после Target_site
  start_idx <- max(1, target_idx - window)
  end_idx <- min(nrow(result), target_idx + window)
  
  result <- result[start_idx:end_idx, , drop = FALSE]
  
  # Сбрасываем rownames
  rownames(result) <- NULL
  
  return(result)
}

# ------------- Function 2: Background threshold:  -------------
getBatchWindows <- function(signals_pos_dye_corrected, batch, window = 50) {
  
  result <- list()

  for (filename in names(signals_pos_dye_corrected)) {
    
    if (signals_pos_dye_corrected[[filename]]$status == "pos_dye") {
      signals <-  signals_pos_dye_corrected[[filename]]$pos_dye_corrected_signals
    } 
    
    if (signals_pos_dye_corrected[[filename]]$status == "run") {
      signals <-  signals_pos_dye_corrected[[filename]]$pos_dye_run_corrected_signals
    }
    
    # Соответствующий read
    read <- batch[[filename]]
    
    # Получаем окно вокруг target site
    result[[filename]] <- getTargetWindow(
      signals_pos_dye_corrected = signals,
      read = read,
      window = window
    )
  }
  
  return(result)
}

# ------------- Function 3: Background threshold: area calculation -------------
trapAreaCalc <- function(sig, peak_start, peak_end) {
  
  pos <- seq(from = peak_start, to = peak_end)
  
  if (length(pos) != length(sig)) {
    stop(paste0("Длина вектора позиций (", length(pos),
                ") не соответствует длине вектора сигналов (", length(sig),")"))
  }
  
  if (length(pos) < 2) {
    return(NA)
  }
  
  dx <- diff(pos)   # расстояние между соседними позициями (убрать нафиг?)
  integral_trap <- sum((sig[-length(sig)] + sig[-1]) / 2 * dx)
  
  return(integral_trap)
}

# ----------- Function 4: Background threshold: run area calculation -----------
areaFracCalculationRun <- function(batch_windows,
                                   batch,
                                   channels = c("A", "C", "G", "T")) {
  
  result <- list()
  
  for (i in seq_along(batch_windows)) {
    
    filename <- batch[[i]]$main_info$filename
    
    window <- batch_windows[[i]]
    n_peaks <- nrow(window)

    areas_summary <- data.frame(
      nuc = window$nuc,
      nt_nuc = window$nt_nuc,
      A_area = numeric(n_peaks),
      C_area = numeric(n_peaks),
      G_area = numeric(n_peaks),
      T_area = numeric(n_peaks)
    ) 
    
    for (j in seq_len(n_peaks)) {
      
      peak_start <- window$peak_start[j]
      peak_end <- window$peak_end[j]
      
      for (channel in channels) {
        sig <- window[[paste0(channel, "_sigs")]][[j]]
        area <- trapAreaCalc(peak_start = peak_start, peak_end = peak_end, sig = sig)
        areas_summary[j, paste0(channel, "_area")] <- area
      }
    }
    
    areas_summary$total_area <- 
      areas_summary$A_area + areas_summary$C_area + 
      areas_summary$G_area + areas_summary$T_area
    
    for (channel in channels) {
      areas_summary[[paste0(channel, "_area_fraction")]] <- 
        areas_summary[[paste0(channel, "_area")]] / areas_summary$total_area
    }
    
    result[[filename]] <- areas_summary
  }
  
  return(result)
}

# ----------- Function 5: Background threshold:  -----------
getMinorValues <- function(area_info, target_site) {
  
  # Убираем Target Site
  background <- area_info[area_info$nt_nuc != target_site, , drop = FALSE]
  
  # Собираем minor peak fractions
  minor_values <- unlist(
    lapply(seq_len(nrow(background)), function(i) {
      nuc <- background$nuc[i]
      minor_channels <- setdiff(c("A", "C", "G", "T"), nuc)
      as.numeric(background[i, paste0(minor_channels, "_area_fraction")])})
  )
  
  # Убираем NA / Inf
  minor_values <- minor_values[is.finite(minor_values)]
  
  return(minor_values)
}

# ----------- Function 6: Background threshold:  -----------
iterativeBackground <- function(minor_values, 
                                z_cutoff = 2.576, 
                                max_iter = 100,
                                tol = 0.01,
                                stable_iter = 2) {
  
  background_values <- minor_values
  
  history <- data.frame(
    iteration = integer(),
    n = integer(),
    mean = numeric(),
    sd = numeric(),
    threshold = numeric(),
    threshold_change = numeric()
  )
  
  stable_count <- 0
  
  for (iter in seq_len(max_iter)) {
    mu <- mean(background_values)
    sd_bg <- sd(background_values)
    
    if (is.na(sd_bg) || sd_bg == 0) {
      break
    }

    threshold <- mu + z_cutoff * sd_bg
    
    # Относительное изменение threshold
    if (iter == 1) {
      threshold_change <- NA
    } else {
      previous_threshold <- history$threshold[nrow(history)]
      threshold_change <- abs(threshold - previous_threshold) / previous_threshold
    }
    
    # Сохраняем историю
    history <- rbind(
      history,
      data.frame(
        iteration = iter,
        n = length(background_values),
        mean = mu,
        sd = sd_bg,
        threshold = threshold,
        threshold_change = threshold_change
      )
    )
    
    # Проверка сходимости
    if (!is.na(threshold_change) && threshold_change < tol) {
      stable_count <- stable_count + 1
    } else {
      stable_count <- 0
    }
    if (stable_count >= stable_iter) {
      break
    }
    
    # Sigma clipping
    keep <- background_values <= threshold
    new_background_values <- background_values[keep]
    
    # Если ничего не удалилось
    if (length(new_background_values) == length(background_values)) {
      break
    }
    
    background_values <- new_background_values
  }
  
  # Финальные параметры
  final_mean <- mean(background_values)
  final_sd <- sd(background_values)
  
  final_threshold <- final_mean + z_cutoff * final_sd
  
  list(
    background_values = background_values,
    threshold = final_threshold,
    mean_background = final_mean,
    sd_background = final_sd,
    n_total = length(minor_values),
    n_background = length(background_values),
    n_removed =
      length(minor_values) -
      length(background_values),
    n_iterations = nrow(history),
    history = history
  )
}

# ----------- Function 7: Background threshold:  -----------
backgroundNoise <- function(area_info, target_site) {
  
  # Получаем minor values
  minor_values <- getMinorValues(area_info = area_info, target_site = target_site)
  
  # Проводим расчет background threshold
  history_iter <- iterativeBackground(minor_values = minor_values)
  
  return(history_iter)
}

# ----------- Function 8: Background threshold:  -----------
backgroundNoiseRun <- function(
    control_area_info,
    sample_area_info,
    control_info,
    sample_info,
    run_correction_required
) {
  
  if (run_correction_required) {
    
    all_minor_values <- c(
      unlist(
        lapply(seq_along(control_area_info), function(i) {
          getMinorValues(
            area_info = control_area_info[[i]],
            target_site = control_info[[i]]$aln_info$ts
          )
        }),
        use.names = FALSE
      ),
      unlist(
        lapply(seq_along(sample_area_info), function(i) {
          getMinorValues(
            area_info = sample_area_info[[i]],
            target_site = sample_info[[i]]$aln_info$ts
          )
        }),
        use.names = FALSE
      )
    )
    
    result <- iterativeBackground(minor_values = all_minor_values)
    
  } else {
    
    calculate_batch <- function(batch_area_info, batch_info) {
      
      result <- list()
      
      for (i in seq_along(batch_area_info)) {
        
        area_info <- batch_area_info[[i]]
        target_site <- batch_info[[i]]$aln_info$ts
        filename <- batch_info[[i]]$main_info$filename
        
        result[[filename]] <- backgroundNoise(
          area_info = area_info,
          target_site = target_site
        )
      }
      
      return(result)
    }
    
    result <- list(
      control = calculate_batch(
        batch_area_info = control_area_info,
        batch_info = control_info
      ),
      sample = calculate_batch(
        batch_area_info = sample_area_info,
        batch_info = sample_info
      )
    )
  }
  
  return(result)
}

# ----------- Function 9: Background threshold:  -----------
backgroundThresholdRun <- function(imported_data, corrected_signals) {
  
  # Data preparation
  control_windows <- getBatchWindows(
    signals_pos_dye_corrected = corrected_signals$control,
    batch = imported_data$control_info,
    window = 50
  )
  
  sample_windows <- getBatchWindows(
    signals_pos_dye_corrected = corrected_signals$sample,
    batch = imported_data$sample_info,
    window = 50
  )
  
  # Area & Fraction calculation
  control_area_info <- areaFracCalculationRun(batch_windows = control_windows, 
                                              batch = imported_data$control_info)
  
  sample_area_info <- areaFracCalculationRun(batch_windows = sample_windows, 
                                             batch = imported_data$sample_info)

  # Check status
  run_correction_required <- any(c(
    vapply(
      corrected_signals$control,
      function(x) identical(x$status, "run"),
      logical(1)
    ),
    vapply(
      corrected_signals$sample,
      function(x) identical(x$status, "run"),
      logical(1)
    )
  ))
  
  # Background Noise
  background_noise <- backgroundNoiseRun(
    control_area_info = control_area_info,
    sample_area_info = sample_area_info,
    control_info = imported_data$control_info,
    sample_info = imported_data$sample_info,
    run_correction_required = run_correction_required
  )
  
  result <- list(
    control = control_area_info,
    sample = sample_area_info,
    background_threshold = background_noise,
    run_correction_required = run_correction_required
  )
  
  return(result)
}

# ==============================================================================
# --------------------------- Mutation Frequency Step --------------------------
# ==============================================================================

# ----------- Function 1: Mutation frequency:  -----------
getTargetFractions <- function(read, window, ts_shift) {
  
  bases <- c("A", "C", "G", "T")
  
  target_site <- read$aln_info$start + ts_shift
  target_row <- window[window$nt_nuc == target_site, , drop = FALSE]
  
  if (nrow(target_row) != 1) {
    stop(paste0("Для read ", read$main_info$filename, " target site найден в ", nrow(target_row), " строках."))
  }
  
  fractions <- c(
    A = as.numeric(target_row$A_area_fraction),
    C = as.numeric(target_row$C_area_fraction),
    G = as.numeric(target_row$G_area_fraction),
    T = as.numeric(target_row$T_area_fraction)
  )
  
  if (length(fractions) != length(bases)) {
    stop(paste0("Количество fractions не соответствует количеству bases для read ", read$main_info$filename))
  }
  
  return(fractions)
}

# ----------- Function 2: Mutation frequency:  -----------
getBT <- function(bt_info, 
                  run_correction_required, 
                  type, 
                  read_name) {
  
  if (run_correction_required) {
    return(as.numeric(bt_info$threshold))
  } else if (!run_correction_required) {
    return(as.numeric(bt_info[[type]][[read_name]]$threshold))
  } else {
    stop(paste0("Неизвестное значение run_correction_required: ", run_correction_required))
  }
}

# ----------- Function 3: Mutation frequency:  -----------
extractMutationRead <- function(read, 
                                area_info, 
                                ts_shift, 
                                bt_info, 
                                run_correction_required, 
                                type, 
                                read_name) {
  
  fractions <- getTargetFractions(
    read = read,
    window = area_info,
    ts_shift = ts_shift
  )
  
  target_site <- read$aln_info$start + ts_shift
  
  target_row <- area_info[
    area_info$nt_nuc == target_site,
    ,
    drop = FALSE
  ]
  
  ref_nuc <- as.character(target_row$nuc[1])
  
  threshold <- getBT(bt_info, run_correction_required, type, read_name)
  
  result <- list(
    fractions = fractions,
    ref_nuc = ref_nuc,
    threshold = threshold,
    run_correction_required = run_correction_required
  )
  
  return(result)
}

# ----------- Function 4: Mutation frequency:  -----------
calculateStats <- function(batch_fractions) {
  
  bases <- c("A", "C", "G", "T")
  
  do.call(rbind, lapply(bases, function(nuc) {
    
    fractions <- vapply(
      batch_fractions,
      function(x) as.numeric(x$fractions[[nuc]]),
      numeric(1)
    )
    
    fractions <- fractions[is.finite(fractions)]
    
    n <- length(fractions)
    
    data.frame(
      nuc = nuc,
      mean_fraction = if (n > 0) mean(fractions) else NA_real_,
      sd_fraction = if (n > 1) sd(fractions) else NA_real_,
      n = n
    )
  }))
}

# ----------- Function 5: Mutation frequency:  -----------
classifyRead <- function(read_data) {
  
  bases <- c("A", "C", "G", "T")
  
  fractions <- read_data$fractions
  threshold <- as.numeric(read_data$threshold)
  
  status <- vapply(bases, function(nuc) {
    
    fraction <- as.numeric(fractions[[nuc]])
    
    if (!is.finite(fraction)) {
      return("undefined")
    }
    if (fraction > threshold) {
      if (nuc == read_data$ref_nuc) {
        return("reference nucleotide")
      }
      return("mutation")
    }
    "background noise"
  }, character(1))
  
  names(status) <- bases
  
  return(status)
}

# ----------- Function 6: Mutation frequency:  -----------
mutationFrequency <- function(result) {
  
  result$MF <- NA_real_
  result$MF_error <- NA_real_
  
  for (i in seq_len(nrow(result))) {
    status <- result$status[i]
    # Reference nucleotide: MF is not calculated
    if (status == "reference nucleotide") {
      result$MF[i] <- NA_real_
      result$MF_error[i] <- NA_real_
      next
    }
    # Background noise: MF is not calculated
    if (status == "background noise") {
      next
    }
    # Undefined status: MF is not calculated
    if (status == "undefined") {
      next
    }
    # Calculate MF only for mutations
    if (status == "mutation") {
      sample_mean <- result$mean_sample[i]
      control_mean <- result$mean_control[i]
      sample_sd <- result$sd_sample[i]
      control_sd <- result$sd_control[i]
      
      if (!is.finite(sample_mean) || !is.finite(control_mean) || 
          sample_mean <= 0 || control_mean <= 0) {
        next
      }
      
      MF <- sample_mean / control_mean
      result$MF[i] <- MF
      
      # Error cannot be estimated without SD values
      if (!is.finite(sample_sd) || !is.finite(control_sd)) {
        next
      }
      
      # Number of reads is needed for standard errors.
      # If n_control and n_sample are not available in result,
      # SD values alone cannot determine the standard errors.
      if (!all(c("n_control", "n_sample") %in% names(result))) {
        next
      }
      n_sample <- result$n_sample[i]
      n_control <- result$n_control[i]
      if (!is.finite(n_sample) || !is.finite(n_control) || n_sample < 1 || n_control < 1) {
        next
      }
      
      sample_se <- sample_sd / sqrt(n_sample)
      control_se <- control_sd / sqrt(n_control)
      
      result$MF_error[i] <- MF * sqrt((sample_se / sample_mean)^2 + (control_se / control_mean)^2)
    }
  }
  
  # Use character columns to represent "-" for reference nucleotide
  result$MF <- as.character(result$MF)
  result$MF_error <- as.character(result$MF_error)
  
  result$MF[result$status == "reference nucleotide"] <- "-"
  result$MF_error[result$status == "reference nucleotide"] <- "-"
  
  return(result)
}

# ----------- Function 7: Mutation frequency:  -----------
identifyMutationType <- function(mutation_table) {
  
  # Достаём референсный нуклеотид
  reference_row <- mutation_table[
    mutation_table$status == "reference nucleotide",
    ,
    drop = FALSE
  ]
  
  if (nrow(reference_row) != 1) {
    stop("Не удалось однозначно определить референсный нуклеотид.")
  }
  
  reference_nuc <- reference_row$nuc
  
  # Списки переходов и трансверсий
  transitions <- c("A>G", "C>T", "G>A", "T>C")
  transversions <- c("A>C", "A>T", "C>A", "C>G", "G>C", "G>T", "T>A", "T>G")
  
  mutation_rows <- which(mutation_table$status == "mutation")
  
  for (i in mutation_rows) {
    
    mutant_nuc <- mutation_table$nuc[i]
    mutation <- paste0(reference_nuc, ">", mutant_nuc)
    
    # Определяем тип
    if (mutation %in% transitions) {
      mutation_table$status[i] <- paste0(mutation, " transition")
      
    } else if (mutation %in% transversions) {
      mutation_table$status[i] <- paste0(mutation,
                                         " transversion")
    }
  }
  
  return(mutation_table)
}

# ----------- Function 8: Mutation frequency:  -----------
mutationFrequencyRun <- function(imported_data, background_threshold, ts_shift) {
  
  bases <- c("A", "C", "G", "T")
  
  control_info <- imported_data$control_info
  sample_info <- imported_data$sample_info
  
  control_area_info <- background_threshold$control
  sample_area_info <- background_threshold$sample
  
  bt_info <- background_threshold$background_threshold
  run_correction_required <- background_threshold$run_correction_required
  
  # 1. Extract fractions and reference nucleotides without classification
  processBatch <- function(batch_info, batch_area_info, type) {
    
    lapply(names(batch_info), function(filename) {
      
      extractMutationRead(
        read = batch_info[[filename]],
        area_info = batch_area_info[[filename]],
        ts_shift = ts_shift,
        bt_info = bt_info,
        run_correction_required = run_correction_required,
        type = type,
        read_name = filename
      )
      
    }) |>
      setNames(names(batch_info))
  }
  
  control_fractions <- processBatch(
    batch_info = control_info,
    batch_area_info = control_area_info,
    type = "control"
  )
  
  sample_fractions <- processBatch(
    batch_info = sample_info,
    batch_area_info = sample_area_info,
    type = "sample"
  )
  
  # 2. Calculate mean fractions + SD
  control_stats <- calculateStats(control_fractions)
  sample_stats <- calculateStats(sample_fractions)
  
  result <- data.frame(
    nuc = bases,
    mean_control = control_stats$mean_fraction,
    sd_control = control_stats$sd_fraction,
    n_control = control_stats$n,
    mean_sample = sample_stats$mean_fraction,
    sd_sample = sample_stats$sd_fraction,
    n_sample = sample_stats$n,
    status = NA_character_,
    stringsAsFactors = FALSE
  )
  
  # 3. Classification
  all_reads <- c(control_fractions, sample_fractions)
  thresholds <- vapply(
    all_reads,
    function(x) as.numeric(x$threshold),
    numeric(1)
  )
  
  ref_nucs <- unique(vapply(
    all_reads,
    function(x) as.character(x$ref_nuc),
    character(1)
  ))
  
  ref_nuc <- ref_nucs
  
  if (!run_correction_required) {
    control_status <- lapply(control_fractions, classifyRead)
    sample_status <- lapply(sample_fractions, classifyRead)
    
    for (i in seq_along(bases)) {
      nuc <- bases[i]
      statuses <- vapply(
        sample_status,
        function(x) x[[nuc]],
        character(1)
      )
      if (any(statuses == "reference nucleotide")) {
        result$status[i] <- "reference nucleotide"
      } else if (any(statuses == "mutation")) {
        result$status[i] <- "mutation"
      } else if (all(statuses == "background noise")) {
        result$status[i] <- "background noise"
      } else {
        result$status[i] <- "undefined"
      }
    }
  } else {
    thresholds <- thresholds[is.finite(thresholds)]
    BT <- unique(thresholds)
    
    for (i in seq_along(bases)) {
      
      nuc <- bases[i]
      
      mean_fraction <- result$mean_sample[i]
      sd_fraction <- result$sd_sample[i]
      
      if (!is.finite(mean_fraction)) {
        result$status[i] <- "undefined"
        next
      }
      
      # If there is only one read, SD is NA.
      lower_bound <- if (is.finite(sd_fraction)) {
        mean_fraction - sd_fraction
      } else {
        mean_fraction
      }
      if (lower_bound > BT) {
        result$status[i] <- if (nuc == ref_nuc) {
          "reference nucleotide"
        } else {
          "mutation"
        }
      } else {
        result$status[i] <- "background noise"
      }
    }
  }
  
  # 4. Calculate mutation frequencies
  mutation_frequency <- mutationFrequency(result = result)
  
  # 5. Identify mutation types
  identified_mutations <- identifyMutationType(
    mutation_table = mutation_frequency
  )
  
  return(identified_mutations)
}

# ----------- Function 9: Mutation frequency:  -----------
formatMutationFrequency <- function(mutation_results, digits = 4) {
  
  format_value_error <- function(value, error) {
    
    result <- rep(NA_character_, length(value))
    
    valid_value <- !is.na(value)
    
    # Если MF задан, сохраняем его значение
    result[valid_value] <- formatC(
      value[valid_value],
      format = "f",
      digits = digits
    )
    
    # Добавляем ± error только если ошибка не NA
    has_error <- valid_value & !is.na(error)
    
    result[has_error] <- paste0(
      formatC(
        value[has_error],
        format = "f",
        digits = digits
      ),
      " ± ",
      formatC(
        error[has_error],
        format = "f",
        digits = digits
      )
    )
    
    result
  }
  
  format_table <- function(df) {
    
    data.frame(
      Nucleotide = df$nuc,
      
      `Control Fraction` = format_value_error(
        df$mean_control,
        df$sd_control
      ),
      
      `Sample Fraction` = format_value_error(
        df$mean_sample,
        df$sd_sample 
      ),
      
      `Mutation Frequency` = format_value_error(
        suppressWarnings(as.numeric(as.character(df$MF))),
        suppressWarnings(as.numeric(as.character(df$MF_error)))
      ),
      
      `Mutation Status` = df$status,
      
      check.names = FALSE,
      stringsAsFactors = FALSE
    ) |>
      (\(out) {
        
        # Сохраняем нечисловые значения MF, например "-"
        special <- is.na(suppressWarnings(
          as.numeric(as.character(df$MF))
        )) & !is.na(df$MF)
        
        out$`Mutation Frequency`[special] <-
          as.character(df$MF[special])
        
        out
        
      })()
  }
  
  lapply(mutation_results, format_table)
}

# ----------- Function 10: Mutation frequency:  -----------
mutationFrequencyRegion <- function(imported_data, background_threshold) {
  
  # Получаем длину target region
  target_region <- strsplit(as.character(imported_data$reference_target_info$target_info$target_seq), "")[[1]]
  
  result <- list()
  
  for (i in seq_along(target_region)) {
    
    nuc <- target_region[i]
    result[[i]] <- mutationFrequencyRun(
      imported_data = imported_data,
      background_threshold = background_threshold,
      ts_shift = i - 1
    )
  }
  
  names(result) <- target_region
  
  formatted_results <- formatMutationFrequency(result)

  return(formatted_results)
}

################################################################################
#                          REPORT PREPARATION STEP                             #
################################################################################

# ----------- Function 1: Report preparation:  -----------
plotTargetRegion <- function(
    corrected_signals,
    imported_data,
    background_threshold,
    filename = NULL,
    width = 14,
    height = 6,
    dpi = 300
) {
  
  bases <- c("A", "C", "G", "T")
  
  nuc_colors <- c(A = "green4", C = "blue2", G = "orange2", T = "red3")
  
  stage_info <- list(
    raw = list(
      field = "raw_signals",
      label = "Raw",
      linetype = "solid",
      color = "grey50"
    ),
    pos = list(
      field = "pos_corrected_signals",
      label = "SPCF",
      linetype = "dashed"
    ),
    pos_dye = list(
      field = "pos_dye_corrected_signals",
      label = "SPCF + SDCF",
      linetype = "dotdash"
    ),
    pos_dye_run = list(
      field = "pos_dye_run_corrected_signals",
      label = "SPCF + SDCF + SRCF",
      linetype = "dotted"
    )
  )
  
  # Получение последовательности целевого региона
  target_seq <- as.character(imported_data$reference_target_info$target_info$target_seq)
  target_seq <- strsplit(target_seq, "")[[1]]
  
  output <- list()
  
  for (type in c("control", "sample")) {
    
    output_type <- if (type == "control") "Control" else "Sample"
    output[[output_type]] <- list()
    
    reads <- names(corrected_signals[[type]])
    
    for (read_name in reads) {
      
      read_corrected <- corrected_signals[[type]][[read_name]]
      
      # Информация о прочтении
      info_name <- paste0(type, "_info")
      read_info <- imported_data[[info_name]][[read_name]]
      
      if (is.null(read_info) || is.null(read_info$aln_info)) {
        warning("Нет aln_info для ", type, ": ", read_name)
        next
      }
      
      start_pos <- read_info$aln_info$start
      end_pos <- read_info$aln_info$end
      
      # Выбираем строки целевого региона по nt_nuc
      raw_table <- read_corrected$raw_signals
      
      if (is.null(raw_table) || nrow(raw_table) == 0) {
        warning("Нет raw_signals для ", read_name)
        next
      }
      
      region_rows <- which(raw_table$nt_nuc >= start_pos & raw_table$nt_nuc <= end_pos)
      
      if (length(region_rows) == 0) {
        warning("Целевой регион не найден для ", read_name)
        next
      }
      
      region_raw <- raw_table[region_rows, , drop = FALSE]
      
      # Получаем последовательность именно для выбранного региона
      region_seq <- target_seq[seq_len(nrow(region_raw))]
      
      if (length(region_seq) != nrow(region_raw)) {
        warning("Длина последовательности не совпадает с регионом: ", read_name)
        next
      }
      
      # Извлекаем порог фонового шума
      if (background_threshold$run_correction_required) {
        threshold <- background_threshold$background_threshold$threshold
      } else {
        threshold <- background_threshold$background_threshold[[type]][[read_name]]$threshold
      }
      
      
      # 1. Секвенограмма: сбор сигналов всех каналов и этапов
      signal_data <- list()
      k <- 0
      
      for (stage in names(stage_info)) {
        
        field <- stage_info[[stage]]$field
        stage_table <- read_corrected[[field]]
        
        # Пропускаем отсутствующую run-коррекцию
        if (is.null(stage_table)) next
        
        selected <- stage_table[stage_table$nt_nuc >= start_pos & stage_table$nt_nuc <= end_pos, , drop = FALSE]
        
        if (nrow(selected) == 0) next
        
        for (i in seq_len(nrow(selected))) {
          
          row <- selected[i, , drop = FALSE]
          nuc_position <- row$nt_nuc
          
          for (base in bases) {
            
            signal_col <- paste0(base, "_sigs")
            if (!signal_col %in% names(selected)) next
            
            sig <- as.numeric(
              unlist(row[[signal_col]], use.names = FALSE)
            )
            
            if (length(sig) == 0) next
            
            # Реальные координаты сигнала
            x <- seq(from = row$peak_start, length.out = length(sig))
            k <- k + 1
            
            signal_data[[k]] <- data.frame(
              x = x,
              signal = sig,
              channel = base,
              stage = stage,
              nt_nuc = nuc_position
            )
          }
        }
      }
      
      if (length(signal_data) == 0) {
        warning("Не удалось собрать сигналы для ", read_name)
        next
      }
      
      signal_df <- do.call(rbind, signal_data)
      
      # Границы пиков: первая позиция — обе границы; остальные — только peak_end.
      boundary_data <- data.frame(
        x = c(region_raw$peak_start[1], region_raw$peak_end),
        boundary_type = c("start", rep("end", nrow(region_raw)))
      )
      
      # Центры пиков и подписи последовательности
      center_data <- data.frame(
        x = region_raw$real_center_pos,
        label = region_seq
      )
      
      # Позиция целевого сайта
      target_site <- imported_data[[paste0(type, "_info")]][[read_name]]$
        aln_info$ts
      
      # Таблица исходных сигналов
      raw_table <- corrected_signals[[type]][[read_name]]$raw_signals
      
      # Находим пик, соответствующий целевому сайту
      target_row <- raw_table[
        raw_table$nt_nuc == target_site,
      ]
      
      if (nrow(target_row) == 1) {
        target_left <- target_row$peak_end
        target_right <- target_row$peak_start
      } else {
        target_left <- NA_real_
        target_right <- NA_real_
        warning(
          "Не удалось однозначно найти целевой сайт: ",
          read_name
        )
      }
      
      p_signal <- ggplot(
        signal_df,
        aes(
          x = x,
          y = signal,
          color = channel,
          linetype = stage,
          group = interaction(channel, stage, nt_nuc)
        )
      ) +
        geom_area(
          aes(
            fill = channel,
            group = interaction(channel, stage, nt_nuc)
          ),
          alpha = 0.15,
          position = "identity"
        ) +
        geom_line(
          linewidth = 0.35,
          na.rm = TRUE
        ) +
        geom_vline(
          data = boundary_data,
          aes(xintercept = x),
          color = "grey55",
          linewidth = 0.35
        ) +
        geom_vline(
          data = center_data,
          aes(xintercept = x),
          color = "black",
          linewidth = 0.25,
          linetype = "dotted"
        ) +
        scale_color_manual(
          name = "Nucleotide",
          values = nuc_colors,
          breaks = bases
        ) +
        scale_fill_manual(
          name = "Nucleotide",
          values = nuc_colors,
          breaks = bases
        ) +
        scale_linetype_manual(
          name = "Correction",
          values = vapply(
            stage_info,
            function(x) x$linetype,
            character(1)
          ),
          breaks = names(stage_info),
          labels = vapply(
            stage_info,
            function(x) x$label,
            character(1)
          )
        ) +
        guides(
          color = guide_legend(order = 1),
          fill = "none",
          linetype = guide_legend(order = 2)
        ) +
        scale_x_continuous(
          breaks = center_data$x,
          labels = center_data$label,
          expand = c(0, 0)
        ) +
        labs(
          title = "Chromatogram",
          x = "Nucleotide",
          y = "Signal intensity"
        ) +
        theme_bw(base_size = 10) +
        theme(
          text = element_text(size = 12),
          axis.text = element_text(size = 12),
          axis.title = element_text(size = 12),
          plot.title = element_text(size = 12),
          legend.text = element_text(size = 12),
          legend.title = element_text(size = 12),
          axis.text.x = element_text(
            size = 12,
            angle = 0,
            vjust = 0.5,
            hjust = 0.5
          ),
          legend.position = "right"
        )
      
      if (is.finite(target_left) && is.finite(target_right)) {
        p_signal <- p_signal +
          ggplot2::annotate(
            "rect",
            xmin = target_left,
            xmax = target_right - 1,
            ymin = -Inf,
            ymax = Inf,
            fill = NA,
            color = "black",
            linewidth = 0.8,
            linetype = "solid"
          )
      }
      
      # 2. Boxplot долей площадей по нуклеотидам
      area_table <- background_threshold[[type]][[read_name]]
      
      if (is.null(area_table)) {
        warning("Нет area_info для ", read_name)
        next
      }
      
      area_region <- area_table[
        area_table$nt_nuc >= start_pos &
          area_table$nt_nuc <= end_pos,
      ]
      
      # Преобразуем таблицу в длинный формат
      area_df <- tidyr::pivot_longer(
        area_region,
        cols = c(
          A_area_fraction,
          C_area_fraction,
          G_area_fraction,
          T_area_fraction
        ),
        names_to = "channel",
        values_to = "area_fraction"
      )
      
      # Названия каналов
      area_df$channel <- sub(
        "_area_fraction",
        "",
        area_df$channel
      )
      
      area_df$channel <- factor(
        area_df$channel,
        levels = c("A", "C", "G", "T")
      )
      
      # Порядок позиций
      position_order <- area_region$nt_nuc
      
      # Подписи — только буквы нуклеотидов
      position_labels <- setNames(
        region_seq[
          match(position_order, area_region$nt_nuc)
        ],
        as.character(position_order)
      )
      
      area_df$nt_nuc <- factor(
        area_df$nt_nuc,
        levels = position_order
      )
      
      area_df <- area_df |>
        dplyr::group_by(nt_nuc) |>
        dplyr::arrange(area_fraction, .by_group = TRUE) |>
        dplyr::mutate(
          stack_order = dplyr::row_number()
        ) |>
        dplyr::ungroup()
      
      # Накопительная столбчатая диаграмма
      p_area <- ggplot2::ggplot(
        data = area_df,
        mapping = ggplot2::aes(
          x = nt_nuc,
          y = area_fraction,
          fill = channel,
          group = interaction(nt_nuc, stack_order)
        )
      ) +
        ggplot2::geom_col(
          width = 0.85,
          position = ggplot2::position_stack(reverse = TRUE),
          alpha = 0.5
        ) +
        geom_hline(
          yintercept = threshold,
          color = "black",
          linetype = "longdash",
          linewidth = 0.5
        ) +
        ggplot2::scale_x_discrete(
          labels = position_labels
        ) +
        ggplot2::scale_fill_manual(
          values = c(
            A = "green4",
            C = "blue2",
            G = "orange2",
            T = "red3"
          ),
          breaks = c("A", "C", "G", "T"),
          name = "Nucleotide"
        ) +
        ggplot2::labs(
          x = "Target Region",
          y = "Area Fraction",
          title = "Fraction Distribution"
        ) +
        ggplot2::annotate(
          "text",
          x = 1,
          y = threshold,
          label = paste0("Threshold = ", round(threshold, 3)),
          hjust = 0,
          vjust = -0.5,
          color = "black"
        ) +
        ggplot2::theme_minimal() +
        ggplot2::theme(
          text = ggplot2::element_text(size = 12),
          axis.text = ggplot2::element_text(size = 12),
          axis.title = ggplot2::element_text(size = 12),
          plot.title = ggplot2::element_text(size = 12),
          legend.text = ggplot2::element_text(size = 12),
          legend.title = ggplot2::element_text(size = 12),
          axis.text.x = ggplot2::element_text(
            size = 12,
            angle = 0,
            vjust = 0.5,
            hjust = 0.5
          ),
          legend.position = "none"
        )
      
      # Объединяем два субплота вертикально
      combined_plot <- p_signal / p_area #+
        #plot_layout(heights = c(2, 1))
      
      output[[output_type]][[read_name]] <- list(
        plot = combined_plot
      )
      
      # Сохранение графика при необходимости
      if (!is.null(filename)) {
        out_file <- file.path(
          filename,
          paste0(output_type, "_", read_name, ".png")
        )
        
        if (!dir.exists(filename)) {
          dir.create(filename, recursive = TRUE)
        }
        
        ggplot2::ggsave(
          filename = out_file,
          plot = combined_plot,
          width = width,
          height = height,
          dpi = dpi,
          device = "png"
        )
        
        message("График сохранён: ", normalizePath(out_file))
      }
    }
  }
  
  return(output)
}


# ----------- Function 2: Report preparation:  -----------
reportPreparationRun <- function(
  imported_data, # from step 3  
  corrected_signals, # from step 4
  background_threshold, # from step 5
  mutation_frequency # from step 6
) {
  
  plots <- plotTargetRegion(
    corrected_signals = corrected_signals,
    imported_data = imported_data,
    background_threshold = background_threshold,
    filename = NULL
    )
  
  return(plots)
}








