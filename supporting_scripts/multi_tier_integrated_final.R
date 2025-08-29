# Multi-Tier Assignment with Signal Balancing Integration (FINAL)
# 
# Author: Manus AI
# Version: Integrated Final 3.0
# Purpose: Combines optimized 5-tier system with signal balancing for strong/weak signal cell types
# 
# NEW FEATURES:
# 1. Strong signal suppression (e.g., Somite over-assignment)
# 2. Weak signal protection (e.g., Gut_surface_ectoderm under-assignment)
# 3. Flexible proportion constraints per cell type
# 4. Enhanced spatial consensus with signal-specific boosts
# 5. All optimizations from previous version maintained

# Required libraries
library(Matrix)

#' Multi-Tier Assignment with Integrated Signal Balancing (FINAL)
#' 
#' Combines the optimized 5-tier system with signal balancing capabilities
#' Solves problems like Somite over-assignment and weak signal under-assignment
#' 
#' @param weights_matrix Normalized RCTD weights matrix (spots x cell_types)
#' @param cell_type_names Vector of cell type names corresponding to matrix columns
#' @param spatial_coords Optional spatial coordinates matrix (spots x 2)
#' 
#' # Original 5-tier parameters (fully preserved)
#' @param excluded_cell_types Vector of cell type names to exclude from assignment
#' @param protected_cell_types Vector of cell type names to protect in Tier 0
#' @param max_proportion Maximum allowed proportion for any single cell type (default: 0.6)
#' @param min_proportion Minimum required proportion for cell types (default: 0.01)
#' @param proportion_correction Whether to apply global proportion correction (default: TRUE)
#' @param correction_method Method for proportion correction (default: "iterative")
#' @param k_neighbors Number of neighbors for spatial consensus (default: 5)
#' @param spatial_weight Weight for spatial information (default: 0.6)
#' @param protection_threshold Minimum weight to protect specified cell types (default: 0.2)
#' @param dominance_threshold Minimum dominance ratio to protect (default: 1.2)
#' @param enable_tier0_protection Whether to enable Tier 0 protection (default: TRUE)
#' 
#' # NEW: Signal balancing parameters
#' @param suppressed_cell_types Vector of cell types to suppress (e.g., c("Somite"))
#' @param max_suppressed_proportion Maximum proportion for suppressed cell types (default: 0.25)
#' @param suppression_strength Factor to reduce suppressed cell type weights (default: 0.8)
#' @param weak_signal_protected_types Vector of weak signal cell types to protect
#' @param weak_signal_threshold Lower protection threshold for weak signals (default: 0.08)
#' @param weak_signal_dominance_threshold Lower dominance requirement for weak signals (default: 1.05)
#' @param weak_signal_spatial_boost Spatial weight boost for weak signals (default: 1.2)
#' @param enable_signal_balancing Whether to enable signal balancing features (default: TRUE)
#' 
#' # Performance parameters (preserved)
#' @param enable_optimizations Whether to use computational optimizations (default: TRUE)
#' @param verbose Whether to print progress messages (default: TRUE)
#' 
#' @return Data frame with assignment results
#' 
#' @export
multi_tier_assignment_integrated_final <- function(weights_matrix,
                                                   cell_type_names,
                                                   spatial_coords = NULL,
                                                   
                                                   # Original 5-tier parameters (fully preserved)
                                                   excluded_cell_types = NULL,
                                                   protected_cell_types = NULL,
                                                   max_proportion = 0.6,
                                                   min_proportion = 0.01,
                                                   proportion_correction = TRUE,
                                                   correction_method = "iterative",
                                                   k_neighbors = 5,
                                                   spatial_weight = 0.6,
                                                   protection_threshold = 0.2,
                                                   dominance_threshold = 1.2,
                                                   enable_tier0_protection = TRUE,
                                                   
                                                   # NEW: Signal balancing parameters
                                                   suppressed_cell_types = NULL,
                                                   max_suppressed_proportion = 0.25,
                                                   suppression_strength = 0.8,
                                                   weak_signal_protected_types = NULL,
                                                   weak_signal_threshold = 0.08,
                                                   weak_signal_dominance_threshold = 1.05,
                                                   weak_signal_spatial_boost = 1.2,
                                                   enable_signal_balancing = TRUE,
                                                   
                                                   # Performance parameters (preserved)
                                                   enable_optimizations = TRUE,
                                                   verbose = TRUE) {
  
  start_time <- Sys.time()
  
  if (verbose) {
    cat("=== Multi-Tier Assignment with Integrated Signal Balancing (FINAL) ===\n")
    cat("Optimizations:", ifelse(enable_optimizations, "ENABLED", "DISABLED"), "\n")
    
    if (enable_tier0_protection && !is.null(protected_cell_types)) {
      cat("Tier 0 protection: ENABLED for", paste(protected_cell_types, collapse = ", "), "\n")
    }
    
    if (enable_signal_balancing) {
      cat("Signal balancing: ENABLED\n")
      if (!is.null(suppressed_cell_types)) {
        cat("Suppressed cell types:", paste(suppressed_cell_types, collapse = ", "), 
            "- max", max_suppressed_proportion * 100, "%\n")
      }
      if (!is.null(weak_signal_protected_types)) {
        cat("Weak signal protected:", paste(weak_signal_protected_types, collapse = ", "), "\n")
      }
    } else {
      cat("Signal balancing: DISABLED\n")
    }
    cat("\n")
  }
  
  # OPTIMIZATION 1: Fast input validation and conversion
  if (enable_optimizations) {
    weights_matrix <- as.matrix(weights_matrix)
    if (!is.null(spatial_coords)) spatial_coords <- as.matrix(spatial_coords)
    
    # Fast NA/Inf handling
    weights_matrix[!is.finite(weights_matrix)] <- 0
    
    # Fast normalization using rowSums
    row_sums <- rowSums(weights_matrix)
    row_sums[row_sums == 0] <- 1
    weights_matrix <- weights_matrix / row_sums
  } else {
    # Original slower method
    if (!is.matrix(weights_matrix)) weights_matrix <- as.matrix(weights_matrix)
    weights_matrix[is.na(weights_matrix)] <- 0
    weights_matrix[is.infinite(weights_matrix)] <- 0
    row_sums <- rowSums(weights_matrix)
    row_sums[row_sums == 0] <- 1
    weights_matrix <- weights_matrix / row_sums
  }
  
  n_spots <- nrow(weights_matrix)
  n_cell_types <- ncol(weights_matrix)
  
  if (verbose) {
    cat("Processing", n_spots, "spots with", n_cell_types, "cell types\n")
  }
  
  # Validate cell type names
  if (length(cell_type_names) != n_cell_types) {
    stop("Length of cell_type_names must match number of columns in weights_matrix")
  }
  colnames(weights_matrix) <- cell_type_names
  
  # STEP 1: Apply exclusion constraints (original logic preserved)
  constrained_weights <- weights_matrix
  constraint_applied <- rep(FALSE, n_spots)
  
  if (!is.null(excluded_cell_types)) {
    excluded_indices <- which(cell_type_names %in% excluded_cell_types)
    
    if (length(excluded_indices) > 0) {
      if (verbose) {
        if (enable_optimizations) {
          # Fast max calculation
          max_indices <- max.col(weights_matrix, ties.method = "first")
          excluded_count <- sum(max_indices %in% excluded_indices)
        } else {
          max_indices <- apply(weights_matrix, 1, which.max)
          excluded_count <- sum(max_indices %in% excluded_indices)
        }
        cat("Excluding", paste(excluded_cell_types, collapse = ", "), 
            "- affects", excluded_count, "spots\n")
      }
      
      # Vectorized exclusion
      constrained_weights[, excluded_indices] <- 0
      
      # Fast renormalization
      if (enable_optimizations) {
        remaining_sums <- rowSums(constrained_weights)
        valid_spots <- remaining_sums > 0
        constrained_weights[valid_spots, ] <- constrained_weights[valid_spots, ] / remaining_sums[valid_spots]
        constraint_applied[valid_spots] <- TRUE
      } else {
        remaining_sums <- rowSums(constrained_weights)
        valid_spots <- remaining_sums > 0
        if (sum(valid_spots) > 0) {
          constrained_weights[valid_spots, ] <- constrained_weights[valid_spots, ] / remaining_sums[valid_spots]
          constraint_applied[valid_spots] <- TRUE
        }
      }
    }
  }
  
  # STEP 2: NEW - Apply signal balancing constraints
  signal_balancing_applied <- rep(FALSE, n_spots)
  
  if (enable_signal_balancing) {
    # Apply strong signal suppression
    if (!is.null(suppressed_cell_types)) {
      if (verbose) cat("Applying strong signal suppression...\n")
      
      for (suppressed_type in suppressed_cell_types) {
        if (suppressed_type %in% cell_type_names) {
          suppressed_idx <- which(cell_type_names == suppressed_type)
          
          # Count spots where this type would be max before suppression
          original_max_indices <- max.col(constrained_weights, ties.method = "first")
          affected_count <- sum(original_max_indices == suppressed_idx)
          
          # Apply suppression
          constrained_weights[, suppressed_idx] <- constrained_weights[, suppressed_idx] * suppression_strength
          signal_balancing_applied <- signal_balancing_applied | (original_max_indices == suppressed_idx)
          
          if (verbose) {
            cat("Suppressed", suppressed_type, "by factor", suppression_strength, 
                "- affects", affected_count, "spots\n")
          }
        }
      }
      
      # Renormalize after suppression
      if (enable_optimizations) {
        remaining_sums <- rowSums(constrained_weights)
        valid_spots <- remaining_sums > 0
        constrained_weights[valid_spots, ] <- constrained_weights[valid_spots, ] / remaining_sums[valid_spots]
      } else {
        remaining_sums <- rowSums(constrained_weights)
        valid_spots <- remaining_sums > 0
        if (sum(valid_spots) > 0) {
          constrained_weights[valid_spots, ] <- constrained_weights[valid_spots, ] / remaining_sums[valid_spots]
        }
      }
    }
  }
  
  # OPTIMIZATION 2: Pre-compute all weight statistics (preserved)
  if (enable_optimizations) {
    # Vectorized max weights calculation
    max_weights <- apply(constrained_weights, 1, max)
    max_indices <- max.col(constrained_weights, ties.method = "first")
    
    # Vectorized second max calculation
    second_max_weights <- numeric(n_spots)
    for (i in 1:n_spots) {
      spot_weights <- constrained_weights[i, ]
      spot_weights[max_indices[i]] <- 0
      second_max_weights[i] <- max(spot_weights)
    }
    
    # Pre-compute weight gaps
    weight_gaps <- max_weights - second_max_weights
  } else {
    # Original method (preserved)
    max_weights <- apply(constrained_weights, 1, function(x) {
      x[is.na(x)] <- 0
      if (all(x == 0)) return(0)
      return(max(x))
    })
    
    max_indices <- apply(constrained_weights, 1, function(x) {
      x[is.na(x)] <- 0
      if (all(x == 0)) return(1)
      return(which.max(x))
    })
    
    second_max_weights <- numeric(n_spots)
    for (i in 1:n_spots) {
      spot_weights <- constrained_weights[i, ]
      spot_weights[is.na(spot_weights)] <- 0
      spot_weights[max_indices[i]] <- 0
      second_max_weights[i] <- max(spot_weights, na.rm = TRUE)
    }
    second_max_weights[is.na(second_max_weights)] <- 0
    
    weight_gaps <- max_weights - second_max_weights
    weight_gaps[is.na(weight_gaps)] <- 0
  }
  
  # Initialize results (preserved structure with new fields)
  spot_ids <- if (is.null(rownames(weights_matrix))) {
    paste0("Spot_", 1:n_spots)
  } else {
    rownames(weights_matrix)
  }
  
  results <- data.frame(
    spot_id = spot_ids,
    original_assignment = character(n_spots),
    assignment_method = character(n_spots),
    assignment_confidence = numeric(n_spots),
    assignment_tier = character(n_spots),
    tier0_protected = rep(FALSE, n_spots),
    protected_cell_type = character(n_spots),
    # NEW: Signal balancing tracking
    signal_balancing_applied = signal_balancing_applied,
    weak_signal_protected = rep(FALSE, n_spots),
    stringsAsFactors = FALSE
  )
  
  # TIER 0: Enhanced protection with weak signal support
  if (enable_tier0_protection && (!is.null(protected_cell_types) || 
                                  (enable_signal_balancing && !is.null(weak_signal_protected_types)))) {
    if (verbose) cat("Applying Tier 0: Enhanced protection with signal balancing...\n")
    
    # Standard protection (preserved)
    if (!is.null(protected_cell_types)) {
      protected_indices <- which(cell_type_names %in% protected_cell_types)
      valid_protected_types <- cell_type_names[protected_indices]
      
      if (length(protected_indices) > 0) {
        total_protected_count <- 0
        
        for (j in 1:length(protected_indices)) {
          pct_idx <- protected_indices[j]
          pct_name <- valid_protected_types[j]
          
          if (enable_optimizations) {
            # VECTORIZED PROTECTION CRITERIA (preserved)
            pct_weights <- constrained_weights[, pct_idx]
            
            protection_mask <- (
              (pct_weights >= protection_threshold) &
              (max_indices == pct_idx) &
              (pct_weights >= second_max_weights * dominance_threshold)
            )
          } else {
            # Original method (preserved)
            pct_weights <- constrained_weights[, pct_idx]
            
            protection_mask <- (
              (pct_weights >= protection_threshold) &
              (max_indices == pct_idx) &
              (pct_weights >= second_max_weights * dominance_threshold)
            )
          }
          
          protection_mask[is.na(protection_mask)] <- FALSE
          protected_count <- sum(protection_mask)
          
          if (verbose) {
            cat("Protected", pct_name, ":", protected_count, "spots (", 
                round(protected_count / n_spots * 100, 1), "%)\n")
          }
          
          if (protected_count > 0) {
            protected_indices_for_this_type <- which(protection_mask)
            results$original_assignment[protected_indices_for_this_type] <- pct_name
            results$assignment_method[protected_indices_for_this_type] <- "Tier0_Protection"
            results$assignment_confidence[protected_indices_for_this_type] <- pct_weights[protected_indices_for_this_type]
            results$assignment_tier[protected_indices_for_this_type] <- "Tier0_Protected"
            results$tier0_protected[protected_indices_for_this_type] <- TRUE
            results$protected_cell_type[protected_indices_for_this_type] <- pct_name
            
            total_protected_count <- total_protected_count + protected_count
          }
        }
        
        if (verbose && total_protected_count > 0) {
          cat("Tier 0 standard protection:", total_protected_count, "spots\n")
        }
      }
    }
    
    # NEW: Weak signal protection
    if (enable_signal_balancing && !is.null(weak_signal_protected_types)) {
      weak_protected_indices <- which(cell_type_names %in% weak_signal_protected_types)
      valid_weak_types <- cell_type_names[weak_protected_indices]
      
      if (length(weak_protected_indices) > 0) {
        total_weak_protected_count <- 0
        
        for (j in 1:length(weak_protected_indices)) {
          weak_idx <- weak_protected_indices[j]
          weak_name <- valid_weak_types[j]
          
          # More lenient protection criteria for weak signals
          weak_weights <- constrained_weights[, weak_idx]
          
          weak_protection_mask <- (
            (weak_weights >= weak_signal_threshold) &
            (max_indices == weak_idx) &
            (weak_weights >= second_max_weights * weak_signal_dominance_threshold) &
            (results$original_assignment == "")  # Only protect unassigned spots
          )
          
          weak_protection_mask[is.na(weak_protection_mask)] <- FALSE
          weak_protected_count <- sum(weak_protection_mask)
          
          if (verbose) {
            cat("Weak signal protected", weak_name, ":", weak_protected_count, "spots (", 
                round(weak_protected_count / n_spots * 100, 1), "%)\n")
          }
          
          if (weak_protected_count > 0) {
            weak_protected_indices_for_this_type <- which(weak_protection_mask)
            results$original_assignment[weak_protected_indices_for_this_type] <- weak_name
            results$assignment_method[weak_protected_indices_for_this_type] <- "Tier0_WeakSignal"
            results$assignment_confidence[weak_protected_indices_for_this_type] <- weak_weights[weak_protected_indices_for_this_type]
            results$assignment_tier[weak_protected_indices_for_this_type] <- "Tier0_WeakSignal"
            results$tier0_protected[weak_protected_indices_for_this_type] <- TRUE
            results$weak_signal_protected[weak_protected_indices_for_this_type] <- TRUE
            results$protected_cell_type[weak_protected_indices_for_this_type] <- weak_name
            
            total_weak_protected_count <- total_weak_protected_count + weak_protected_count
          }
        }
        
        if (verbose && total_weak_protected_count > 0) {
          cat("Tier 0 weak signal protection:", total_weak_protected_count, "spots\n")
        }
      }
    }
  }
  
  # TIER 1-5: Optimized 5-tier system (preserved with enhancements)
  if (verbose) cat("Applying optimized 5-tier system to remaining spots...\n")
  
  remaining_mask <- results$original_assignment == ""
  remaining_indices <- which(remaining_mask)
  remaining_count <- length(remaining_indices)
  
  if (remaining_count > 0) {
    if (verbose) cat("Processing", remaining_count, "remaining spots\n")
    
    # TIER 1-2: Optimized confidence-based assignment (preserved)
    if (enable_optimizations) {
      # Pre-filter data for remaining spots
      remaining_max_weights <- max_weights[remaining_indices]
      remaining_weight_gaps <- weight_gaps[remaining_indices]
      remaining_max_indices <- max_indices[remaining_indices]
      
      # TIER 1: Vectorized high confidence
      tier1_mask <- (remaining_weight_gaps >= 0.3) | (remaining_max_weights >= 0.5)
      tier1_count <- sum(tier1_mask)
      
      if (tier1_count > 0) {
        tier1_global_indices <- remaining_indices[tier1_mask]
        results$original_assignment[tier1_global_indices] <- cell_type_names[remaining_max_indices[tier1_mask]]
        results$assignment_method[tier1_global_indices] <- "High_Confidence"
        results$assignment_confidence[tier1_global_indices] <- remaining_max_weights[tier1_mask]
        results$assignment_tier[tier1_global_indices] <- "Tier1_High_Confidence"
        
        # Update remaining
        remaining_mask[tier1_global_indices] <- FALSE
        remaining_indices <- which(remaining_mask)
        remaining_count <- length(remaining_indices)
      }
      
      # TIER 2: Vectorized medium confidence
      if (remaining_count > 0) {
        remaining_max_weights <- max_weights[remaining_indices]
        remaining_weight_gaps <- weight_gaps[remaining_indices]
        remaining_max_indices <- max_indices[remaining_indices]
        
        tier2_mask <- (remaining_weight_gaps >= 0.15) | (remaining_max_weights >= 0.4)
        tier2_count <- sum(tier2_mask)
        
        if (tier2_count > 0) {
          tier2_global_indices <- remaining_indices[tier2_mask]
          results$original_assignment[tier2_global_indices] <- cell_type_names[remaining_max_indices[tier2_mask]]
          results$assignment_method[tier2_global_indices] <- "Medium_Confidence"
          results$assignment_confidence[tier2_global_indices] <- remaining_max_weights[tier2_mask]
          results$assignment_tier[tier2_global_indices] <- "Tier2_Medium_Confidence"
          
          # Update remaining
          remaining_mask[tier2_global_indices] <- FALSE
          remaining_indices <- which(remaining_mask)
          remaining_count <- length(remaining_indices)
        }
      }
      
      if (verbose) cat("Tier 1-2: Assigned", tier1_count + tier2_count, "high/medium confidence spots\n")
    }
    
    # TIER 3: Enhanced spatial consensus with weak signal boost
    if (!is.null(spatial_coords) && remaining_count > 0) {
      if (verbose) cat("Tier 3: Applying enhanced spatial consensus with signal balancing...\n")
      
      if (enable_optimizations) {
        spatial_results <- perform_enhanced_spatial_consensus(
          constrained_weights[remaining_indices, , drop = FALSE],
          cell_type_names,
          spatial_coords[remaining_indices, , drop = FALSE],
          spatial_weight,
          k_neighbors,
          weak_signal_protected_types,
          weak_signal_spatial_boost,
          enable_signal_balancing
        )
      } else {
        spatial_results <- perform_spatial_consensus_safe(
          constrained_weights[remaining_indices, , drop = FALSE],
          cell_type_names,
          spatial_coords[remaining_indices, , drop = FALSE],
          spatial_weight,
          k_neighbors
        )
      }
      
      results$original_assignment[remaining_indices] <- spatial_results$assignments
      results$assignment_method[remaining_indices] <- "Spatial_Consensus"
      results$assignment_confidence[remaining_indices] <- spatial_results$confidence
      results$assignment_tier[remaining_indices] <- "Tier3_Spatial_Consensus"
      
      if (verbose) cat("Tier 3: Assigned", remaining_count, "spots using enhanced spatial consensus\n")
      
      # Update remaining
      remaining_mask <- results$original_assignment == ""
      remaining_indices <- which(remaining_mask)
      remaining_count <- length(remaining_indices)
    }
    
    # TIER 4-5: Vectorized weighted and forced assignment (preserved)
    if (remaining_count > 0) {
      if (enable_optimizations) {
        # Vectorized final assignments
        remaining_weights <- constrained_weights[remaining_indices, , drop = FALSE]
        
        # Tier 4: Enhanced weighted assignment
        squared_weights <- remaining_weights^2
        row_sums_sq <- rowSums(squared_weights)
        row_sums_sq[row_sums_sq == 0] <- 1
        squared_weights <- squared_weights / row_sums_sq
        
        tier4_max_indices <- max.col(squared_weights, ties.method = "first")
        tier4_confidence <- squared_weights[cbind(1:remaining_count, tier4_max_indices)]
        
        results$original_assignment[remaining_indices] <- cell_type_names[tier4_max_indices]
        results$assignment_method[remaining_indices] <- "Weighted_Assignment"
        results$assignment_confidence[remaining_indices] <- tier4_confidence
        results$assignment_tier[remaining_indices] <- "Tier4_Weighted"
        
        if (verbose) cat("Tier 4-5: Assigned", remaining_count, "spots using vectorized weighted assignment\n")
      } else {
        # Original method for Tier 4-5 (preserved)
        for (i in 1:remaining_count) {
          global_idx <- remaining_indices[i]
          spot_weights <- constrained_weights[global_idx, ]
          spot_weights[is.na(spot_weights)] <- 0
          
          if (all(spot_weights == 0)) {
            max_idx <- 1
            confidence <- 0
          } else {
            max_idx <- which.max(spot_weights)
            confidence <- spot_weights[max_idx]
          }
          
          results$original_assignment[global_idx] <- cell_type_names[max_idx]
          results$assignment_method[global_idx] <- "Forced_Assignment"
          results$assignment_confidence[global_idx] <- confidence
          results$assignment_tier[global_idx] <- "Tier5_Forced"
        }
        
        if (verbose) cat("Tier 4-5: Assigned", remaining_count, "spots using original method\n")
      }
    }
  }
  
  # Enhanced proportion correction with signal balancing
  if (proportion_correction) {
    if (verbose) cat("Applying enhanced proportion correction with signal balancing...\n")
    
    if (enable_optimizations) {
      corrected_results <- apply_enhanced_proportion_correction(
        results, constrained_weights, cell_type_names, 
        max_proportion, excluded_cell_types, protected_cell_types,
        suppressed_cell_types, max_suppressed_proportion,
        weak_signal_protected_types, enable_signal_balancing, verbose
      )
    } else {
      corrected_results <- apply_proportion_correction_with_tier0_protection(
        results, constrained_weights, cell_type_names, 
        max_proportion, excluded_cell_types, protected_cell_types, correction_method, verbose
      )
    }
    
    final_results <- results
    final_results$final_assignment <- corrected_results$final_assignment
    final_results$corrected <- corrected_results$corrected
    final_results$alternative_assignment <- corrected_results$alternative_assignment
  } else {
    final_results <- results
    final_results$final_assignment <- results$original_assignment
    final_results$corrected <- rep(FALSE, n_spots)
    final_results$alternative_assignment <- rep(NA, n_spots)
  }
  
  # Add weight information (optimized, preserved)
  if (enable_optimizations) {
    for (i in 1:length(cell_type_names)) {
      col_name <- paste0(cell_type_names[i], "_weight")
      final_results[[col_name]] <- weights_matrix[, i]
    }
  } else {
    for (i in 1:length(cell_type_names)) {
      col_name <- paste0(cell_type_names[i], "_weight")
      final_results[[col_name]] <- weights_matrix[, i]
    }
  }
  
  # Add constraint information (enhanced)
  final_results$biological_constraint_applied <- constraint_applied
  final_results$excluded_from_assignment <- sapply(final_results$final_assignment, 
                                                   function(x) x %in% excluded_cell_types)
  
  # Performance reporting (enhanced)
  end_time <- Sys.time()
  elapsed_time <- as.numeric(difftime(end_time, start_time, units = "secs"))
  
  if (verbose) {
    cat("\n=== INTEGRATED FINAL Performance Report ===\n")
    cat("Total processing time:", round(elapsed_time, 2), "seconds\n")
    cat("Processing speed:", round(n_spots / elapsed_time, 0), "spots/second\n")
    cat("Optimizations:", ifelse(enable_optimizations, "ENABLED", "DISABLED"), "\n")
    cat("Signal balancing:", ifelse(enable_signal_balancing, "ENABLED", "DISABLED"), "\n")
    
    generate_integrated_final_report(final_results, excluded_cell_types, protected_cell_types, 
                                    suppressed_cell_types, weak_signal_protected_types,
                                    max_proportion, max_suppressed_proportion, n_spots, 
                                    enable_tier0_protection, enable_signal_balancing)
  }
  
  return(final_results)
}

#' Enhanced spatial consensus with weak signal boost
perform_enhanced_spatial_consensus <- function(weights_matrix, cell_type_names, spatial_coords, 
                                              spatial_weight, k_neighbors,
                                              weak_signal_protected_types = NULL,
                                              weak_signal_spatial_boost = 1.2,
                                              enable_signal_balancing = TRUE) {
  
  weights_matrix <- as.matrix(weights_matrix)
  spatial_coords <- as.matrix(spatial_coords)
  n_spots <- nrow(weights_matrix)
  
  if (n_spots == 0) {
    return(list(assignments = character(0), confidence = numeric(0)))
  }
  
  assignments <- character(n_spots)
  confidence <- numeric(n_spots)
  
  if (n_spots == 1) {
    spot_weights <- weights_matrix[1, ]
    max_idx <- which.max(spot_weights)
    assignments[1] <- cell_type_names[max_idx]
    confidence[1] <- spot_weights[max_idx]
    return(list(assignments = assignments, confidence = confidence))
  }
  
  # OPTIMIZATION: Vectorized distance calculation with signal balancing
  for (i in 1:n_spots) {
    current_coord <- spatial_coords[i, ]
    
    # Vectorized distance calculation
    coord_diff <- sweep(spatial_coords, 2, current_coord, "-")
    distances <- sqrt(rowSums(coord_diff^2))
    distances[i] <- Inf  # Exclude self
    
    # Fast neighbor finding
    neighbor_indices <- order(distances)[1:min(k_neighbors, n_spots-1)]
    
    # Vectorized neighbor weight calculation
    neighbor_weights <- colMeans(weights_matrix[neighbor_indices, , drop = FALSE])
    self_weights <- weights_matrix[i, ]
    
    # Apply signal balancing to spatial weight
    effective_spatial_weight <- spatial_weight
    
    if (enable_signal_balancing && !is.null(weak_signal_protected_types)) {
      # Check if any weak signal type is dominant in neighbors
      max_neighbor_idx <- which.max(neighbor_weights)
      max_neighbor_type <- cell_type_names[max_neighbor_idx]
      
      if (max_neighbor_type %in% weak_signal_protected_types) {
        effective_spatial_weight <- spatial_weight * weak_signal_spatial_boost
        effective_spatial_weight <- min(effective_spatial_weight, 0.95)  # Cap at 95%
      }
    }
    
    # Fast combination with potentially boosted spatial weight
    combined_weights <- (1 - effective_spatial_weight) * self_weights + effective_spatial_weight * neighbor_weights
    
    # Fast assignment
    max_idx <- which.max(combined_weights)
    assignments[i] <- cell_type_names[max_idx]
    confidence[i] <- combined_weights[max_idx]
  }
  
  return(list(assignments = assignments, confidence = confidence))
}

#' Enhanced proportion correction with signal balancing
apply_enhanced_proportion_correction <- function(results, weights_matrix, cell_type_names, 
                                                max_proportion, excluded_cell_types, 
                                                protected_cell_types,
                                                suppressed_cell_types, max_suppressed_proportion,
                                                weak_signal_protected_types, enable_signal_balancing,
                                                verbose) {
  
  n_spots <- nrow(results)
  current_assignments <- results$original_assignment
  
  # Fast proportion calculation
  current_props <- table(current_assignments)
  current_prop_pct <- current_props / n_spots
  
  # Enhanced over-representation detection with signal balancing
  excluded_from_correction <- c(excluded_cell_types, protected_cell_types)
  if (enable_signal_balancing && !is.null(weak_signal_protected_types)) {
    excluded_from_correction <- c(excluded_from_correction, weak_signal_protected_types)
  }
  
  over_represented_types <- c()
  
  # Check general proportion limits
  valid_props <- current_prop_pct[!names(current_prop_pct) %in% excluded_from_correction]
  general_over_represented <- names(valid_props)[valid_props > max_proportion]
  over_represented_types <- c(over_represented_types, general_over_represented)
  
  # Check suppressed cell type limits
  if (enable_signal_balancing && !is.null(suppressed_cell_types)) {
    for (suppressed_type in suppressed_cell_types) {
      if (suppressed_type %in% names(current_prop_pct)) {
        if (current_prop_pct[suppressed_type] > max_suppressed_proportion) {
          over_represented_types <- c(over_represented_types, suppressed_type)
          if (verbose) {
            cat("Suppressed type", suppressed_type, "exceeds limit:", 
                round(current_prop_pct[suppressed_type] * 100, 1), "% >", 
                max_suppressed_proportion * 100, "%\n")
          }
        }
      }
    }
  }
  
  over_represented_types <- unique(over_represented_types)
  
  if (length(over_represented_types) == 0) {
    return(list(
      final_assignment = current_assignments,
      corrected = rep(FALSE, n_spots),
      alternative_assignment = get_optimized_alternative_assignments(weights_matrix, cell_type_names, excluded_cell_types)
    ))
  }
  
  # Enhanced correction with signal balancing awareness
  corrected_assignments <- current_assignments
  corrected_flags <- rep(FALSE, n_spots)
  alternative_assignments <- get_optimized_alternative_assignments(weights_matrix, cell_type_names, excluded_cell_types)
  
  for (over_ct in over_represented_types) {
    # Determine target proportion based on type
    target_proportion <- max_proportion
    if (enable_signal_balancing && !is.null(suppressed_cell_types) && over_ct %in% suppressed_cell_types) {
      target_proportion <- min(max_proportion, max_suppressed_proportion)
    }
    
    current_count <- sum(corrected_assignments == over_ct)
    target_count <- floor(n_spots * target_proportion)
    excess_count <- current_count - target_count
    
    if (excess_count > 0) {
      over_spots <- which(corrected_assignments == over_ct)
      
      # Enhanced protection check
      non_protected_over <- over_spots[
        !results$tier0_protected[over_spots] &
        !results$weak_signal_protected[over_spots] &
        results$assignment_tier[over_spots] != "Tier1_High_Confidence"
      ]
      
      if (length(non_protected_over) > 0) {
        # Vectorized confidence sorting
        spot_confidences <- results$assignment_confidence[non_protected_over]
        spots_to_reassign <- non_protected_over[order(spot_confidences)][1:min(excess_count, length(non_protected_over))]
        
        # Enhanced reassignment with signal balancing awareness
        for (spot_idx in spots_to_reassign) {
          alt_assignment <- alternative_assignments[spot_idx]
          
          # Avoid reassigning to other over-represented types
          if (!alt_assignment %in% over_represented_types) {
            corrected_assignments[spot_idx] <- alt_assignment
            corrected_flags[spot_idx] <- TRUE
          }
        }
        
        if (verbose) {
          actual_reassigned <- sum(corrected_flags[spots_to_reassign])
          cat("Reassigned", actual_reassigned, "spots from", over_ct, "\n")
        }
      }
    }
  }
  
  return(list(
    final_assignment = corrected_assignments,
    corrected = corrected_flags,
    alternative_assignment = alternative_assignments
  ))
}

#' Generate integrated final report
generate_integrated_final_report <- function(final_results, excluded_cell_types, protected_cell_types, 
                                            suppressed_cell_types, weak_signal_protected_types,
                                            max_proportion, max_suppressed_proportion, n_spots, 
                                            enable_tier0_protection, enable_signal_balancing) {
  cat("\n=== Integrated Final Results Summary ===\n")
  cat("Total spots:", n_spots, "\n")
  
  # Tier 0 protection summary
  if (enable_tier0_protection) {
    tier0_count <- sum(final_results$tier0_protected, na.rm = TRUE)
    cat("Tier 0 protected:", tier0_count, "spots\n")
    
    if (!is.null(protected_cell_types)) {
      standard_protected <- sum(final_results$tier0_protected & !final_results$weak_signal_protected, na.rm = TRUE)
      cat("  - Standard protection:", standard_protected, "spots\n")
    }
    
    if (enable_signal_balancing && !is.null(weak_signal_protected_types)) {
      weak_protected <- sum(final_results$weak_signal_protected, na.rm = TRUE)
      cat("  - Weak signal protection:", weak_protected, "spots\n")
    }
  }
  
  # Signal balancing summary
  if (enable_signal_balancing) {
    signal_balanced <- sum(final_results$signal_balancing_applied, na.rm = TRUE)
    cat("Signal balancing applied:", signal_balanced, "spots\n")
  }
  
  # Final proportions with status indicators
  final_props <- table(final_results$final_assignment)
  cat("\nFinal proportions:\n")
  for (i in 1:length(final_props)) {
    cell_type <- names(final_props)[i]
    count <- final_props[i]
    proportion <- count / n_spots
    
    status_indicators <- c()
    if (!is.null(excluded_cell_types) && cell_type %in% excluded_cell_types) {
      status_indicators <- c(status_indicators, "EXCLUDED")
    }
    if (!is.null(protected_cell_types) && cell_type %in% protected_cell_types) {
      status_indicators <- c(status_indicators, "PROTECTED")
    }
    if (enable_signal_balancing && !is.null(suppressed_cell_types) && cell_type %in% suppressed_cell_types) {
      status_indicators <- c(status_indicators, "SUPPRESSED")
    }
    if (enable_signal_balancing && !is.null(weak_signal_protected_types) && cell_type %in% weak_signal_protected_types) {
      status_indicators <- c(status_indicators, "WEAK SIGNAL PROTECTED")
    }
    
    status_str <- if (length(status_indicators) > 0) {
      paste0(" [", paste(status_indicators, collapse = ", "), "]")
    } else {
      ""
    }
    
    cat(cell_type, ":", count, "(", round(proportion * 100, 1), "%)", status_str, "\n")
  }
  
  # Tier distribution
  tier_table <- table(final_results$assignment_tier)
  cat("\nTier distribution:\n")
  print(tier_table)
  
  # Correction summary
  corrected_count <- sum(final_results$corrected, na.rm = TRUE)
  cat("\nSpots corrected by proportion adjustment:", corrected_count, "(", 
      round(corrected_count / n_spots * 100, 1), "%)\n")
}

# Preserve original helper functions for compatibility
get_optimized_alternative_assignments <- function(weights_matrix, cell_type_names, excluded_cell_types) {
  
  n_spots <- nrow(weights_matrix)
  valid_indices <- which(!cell_type_names %in% excluded_cell_types)
  
  if (length(valid_indices) == 0) {
    return(rep(cell_type_names[1], n_spots))
  }
  
  # Vectorized alternative assignment
  valid_weights <- weights_matrix[, valid_indices, drop = FALSE]
  valid_names <- cell_type_names[valid_indices]
  
  # Fast second-best calculation
  alternative_assignments <- character(n_spots)
  
  for (i in 1:n_spots) {
    spot_weights <- valid_weights[i, ]
    sorted_indices <- order(spot_weights, decreasing = TRUE)
    
    if (length(sorted_indices) >= 2) {
      alternative_assignments[i] <- valid_names[sorted_indices[2]]
    } else {
      alternative_assignments[i] <- valid_names[sorted_indices[1]]
    }
  }
  
  return(alternative_assignments)
}

# Fallback functions for compatibility
perform_spatial_consensus_safe <- function(weights_matrix, cell_type_names, spatial_coords, 
                                          spatial_weight, k_neighbors) {
  return(perform_enhanced_spatial_consensus(weights_matrix, cell_type_names, spatial_coords, 
                                           spatial_weight, k_neighbors, NULL, 1.0, FALSE))
}

apply_proportion_correction_with_tier0_protection <- function(results, weights_matrix, cell_type_names, 
                                                             max_proportion, excluded_cell_types, 
                                                             protected_cell_types, correction_method, verbose) {
  return(apply_enhanced_proportion_correction(results, weights_matrix, cell_type_names, 
                                             max_proportion, excluded_cell_types, 
                                             protected_cell_types, NULL, 0.25, NULL, FALSE, verbose))
}

# Quick solution functions for common problems

#' Quick solution for Somite over-assignment problem
solve_somite_overassignment_integrated <- function(weights_matrix,
                                                   cell_type_names,
                                                   spatial_coords,
                                                   somite_max_proportion = 0.25,
                                                   gut_ectoderm_protection = TRUE,
                                                   verbose = TRUE) {
  
  if (verbose) {
    cat("=== Quick Integrated Solution for Somite Over-assignment ===\n")
    cat("Somite max proportion:", somite_max_proportion * 100, "%\n")
    cat("Gut surface ectoderm protection:", gut_ectoderm_protection, "\n\n")
  }
  
  # Set up parameters
  suppressed_types <- c("Somite")
  weak_protected_types <- if (gut_ectoderm_protection) c("Gut_surface_ectoderm") else NULL
  
  # Call integrated function
  results <- multi_tier_assignment_integrated_final(
    weights_matrix = weights_matrix,
    cell_type_names = cell_type_names,
    spatial_coords = spatial_coords,
    
    # Exclude TE (common requirement)
    excluded_cell_types = c("TE"),
    
    # Signal balancing for Somite problem
    suppressed_cell_types = suppressed_types,
    max_suppressed_proportion = somite_max_proportion,
    suppression_strength = 0.8,
    
    # Weak signal protection
    weak_signal_protected_types = weak_protected_types,
    weak_signal_threshold = 0.08,
    weak_signal_dominance_threshold = 1.05,
    weak_signal_spatial_boost = 1.3,
    
    # Enhanced spatial constraints
    spatial_weight = 0.8,
    k_neighbors = 8,
    
    # General parameters
    max_proportion = 0.6,
    enable_signal_balancing = TRUE,
    enable_optimizations = TRUE,
    verbose = verbose
  )
  
  return(results)
}

cat("=== INTEGRATED FINAL Multi-Tier Assignment System Loaded ===\n")
cat("Main function: multi_tier_assignment_integrated_final()\n")
cat("Quick solution: solve_somite_overassignment_integrated()\n")
cat("\nNEW FEATURES:\n")
cat("1. Strong signal suppression (e.g., Somite over-assignment)\n")
cat("2. Weak signal protection (e.g., Gut_surface_ectoderm under-assignment)\n")
cat("3. Flexible proportion constraints per cell type\n")
cat("4. Enhanced spatial consensus with signal-specific boosts\n")
cat("5. All optimizations from previous version maintained\n")
cat("\nUSAGE:\n")
cat("# Basic usage with your existing parameters\n")
cat("results <- multi_tier_assignment_integrated_final(\n")
cat("  weights_matrix, cell_type_names, spatial_coords,\n")
cat("  excluded_cell_types = c('TE'),\n")
cat("  protected_cell_types = c('YS.Endoderm_2'),  # Your existing protection\n")
cat("  suppressed_cell_types = c('Somite'),        # NEW: Suppress Somite\n")
cat("  weak_signal_protected_types = c('Gut_surface_ectoderm')  # NEW: Protect weak signals\n")
cat(")\n")
cat("\n# Quick solution for Somite problem\n")
cat("results <- solve_somite_overassignment_integrated(\n")
cat("  weights_matrix, cell_type_names, spatial_coords,\n")
cat("  somite_max_proportion = 0.25\n")
cat(")\n")

