format_vars <- function(obs_list,
                        agg_levels,
                        poly_static,
                        poly_dynamic,
                        pixel_static,
                        pixel_dynamic,
                        times_dynamic,
                        inc_times,
                        inc_poly,
                        pop_dynamic,
                        pop_static,
                        n_subset=100){
  u_agg_levels <- sort(unique(agg_levels))
  if(length(u_agg_levels) == 1){
    if(u_agg_levels == 2 |u_agg_levels == 3){
      print("all static or all dynamic")
      n <- length(obs_list[[1]])
      index_use <- sample.int(n, n_subset)
      out_list <- lapply(obs_list, function(l) l[index_use])
      pop <- rep(1, n)
      polygon_start_index <- (1:n_subset) - 1
      polygon_end_index <- 1:n_subset
    }
  }else{
    polygon_start_index <- rep(NA, n_subset)
    polygon_end_index <- rep(NA, n_subset)
    if(length(u_agg_levels) == 2){
      # print("test")
      if(all(u_agg_levels == c(2, 3))){
        print("static and dynamic")
        which_static <- which(agg_levels == 2)
        n_static <- length(which_static)
        which_dynamic <- which(agg_levels == 3)
        n_dynamic <- length(which_dynamic)
        
        n <- length(pixel_static)
        index_use <- sample(n, n_subset)
        pixel_use <- pixel_static[index_use]
        index_dynamic <- which(pixel_dynamic %in% pixel_use)
        
        out_list <- list()
        for(i in 1:length(obs_list)){
          out_list[[i]] <- rep(NA, length(index_dynamic))
        }
        
        count <- 0
        for(i in 1:n_subset){
          polygon_start_index[i] <- count
          which_match <- which(pixel_use[i] == pixel_dynamic)
          n_match <- length(which_match)
          which_i <- count + (1:n_match)
          count <- count + n_match
          polygon_end_index[i] <- count
          for(j in 1:n_static){
            out_list[[which_static[j]]][which_i] <- obs_list[[which_static[j]]][index_use[i]]
          }
          for(j in 1:n_dynamic){
            out_list[[which_dynamic[j]]][which_i] <- obs_list[[which_dynamic[j]]][which_match]
          }
        }
        pop <- rep(1, length(index_dynamic))
      }
      
      if(all(u_agg_levels == c(1, 3))){
        print("dynamic and incidence")
        which_inc <- which(agg_levels == 1)
        n_inc <- length(which_inc)
        which_dynamic <- which(agg_levels == 3)
        n_dynamic <- length(which_dynamic)
        
        #choose some polygons
        n_polygons <- length(unique(inc_poly))
        poly_use <- sample(unique(inc_poly), n_subset)
        time_use <- c()
        n_out <- 0
        for(i in 1:n_subset){
          which_times <- unique(inc_times[inc_poly == poly_use[i]])
          time_use[i] <- sample(which_times, 1)
          n_out <- n_out + sum(poly_dynamic == poly_use[i] & times_dynamic == time_use[i])
        }
        
        # time_use <- sample.int(12, n_subset, replace = TRUE)
        index_use <- cbind(poly_use, time_use)
        
        # index_dynamic <- which(poly_dynamic %in% poly_use)
        
        out_list <- list()
        for(i in 1:length(obs_list)){
          out_list[[i]] <- rep(NA, n_out)
        }
        pop <- rep(NA, n_out)
        
        count <- 0
        for(i in 1:n_subset){
          # which_inc_times_match <- which(inc_times == time_use[i])
          # which_inc_poly_match <- which(inc_poly == poly_use[i])
          # which_inc_match <- intersect(which_inc_times_match, which_inc_poly_match)
          which_inc_match <- which(inc_times == time_use[i] & inc_poly == poly_use[i])
          
          
          # which_dynamic_times_match <- which(times_dynamic == time_use[i])
          # which_dynamic_poly_match <- which(poly_dynamic == poly_use[i])
          # which_dynamic_match <- intersect(which_inc_times_match, which_inc_poly_match)
          which_dynamic_match <- which(times_dynamic == time_use[i] & poly_dynamic == poly_use[i])
          
          polygon_start_index[i] <- count
          n_match <- length(which_dynamic_match)
          which_i <- count + (1:n_match)
          
          count <- count + n_match
          polygon_end_index[i] <- count
          # if(i == 1) print(which_match)
          for(j in 1:n_inc){
            out_list[[which_inc[j]]][which_i] <- obs_list[[which_inc[j]]][which_inc_match]
          }
          for(j in 1:n_dynamic){
            out_list[[which_dynamic[j]]][which_i] <- obs_list[[which_dynamic[j]]][which_dynamic_match]
          }
          pop[which_i] <- pop_dynamic[which_dynamic_match]
        }
      }
      
      if(all(u_agg_levels == c(1, 2))){
        print("static and incidence")
        which_inc <- which(agg_levels == 1)
        n_inc <- length(which_inc)
        which_static <- which(agg_levels == 2)
        n_static <- length(which_static)
        
        #choose some polygons
        n_polygons <- length(unique(inc_poly))
        poly_use <- sample(unique(inc_poly), n_subset)
        
        # time_use <- sample.int(12, n_subset, replace = TRUE)
        index_use <- poly_use
        
        ##calculate the number of observations will have
        poly_size <- sapply(poly_use, function(i) sum(poly_static == i))
        poly_times <- sapply(poly_use, function(i) sum(inc_poly == i))
        
        n_out <- sum(poly_size * poly_times)
        
        out_list <- list()
        for(i in 1:length(obs_list)){
          out_list[[i]] <- rep(NA, n_out)
        }
        pop <- rep(NA, n_out)
        
        count <- 0
        for(i in 1:n_subset){
          which_inc_match <- which(inc_poly == poly_use[i])
          which_static_match <- which(poly_static == poly_use[i])
          
          polygon_start_index[i] <- count
          n_match <- poly_size[i] * poly_times[i]
          which_i <- count + (1:n_match)
          count <- count + n_match
          polygon_end_index[i] <- count
          # if(i == 1) print(which_match)
          for(j in 1:n_inc){
            out_list[[which_inc[j]]][which_i] <- rep(obs_list[[which_inc[j]]][which_inc_match],
                                                     poly_size[i])
          }
          for(j in 1:n_static){
            out_list[[which_static[j]]][which_i] <- rep(obs_list[[which_static[j]]][which_static_match],
                                                        each=poly_times[i])
          }
          pop[which_i] <- rep(pop_dynamic[which_static_match],
                              each=poly_times[i])
        }
      }
    }else{
      print("inc, static and dynamic")
      which_inc <- which(agg_levels == 1)
      n_inc <- length(which_inc)
      which_static <- which(agg_levels == 2)
      n_static <- length(which_static)
      which_dynamic <- which(agg_levels == 3)
      n_dynamic <- length(which_dynamic)
      
      #choose some polygons
      n_polygons <- length(unique(inc_poly))
      poly_use <- sample(unique(inc_poly), n_subset)
      
      index_use <- poly_use
      
      ##calculate the number of observations will have
      poly_size <- sapply(poly_use, function(i) sum(poly_static == i))
      poly_times <- sapply(poly_use, function(i) sum(inc_poly == i))
      
      n_out <- sum(poly_size * poly_times)
      
      
      
      out_list <- list()
      for(i in 1:length(obs_list)){
        out_list[[i]] <- rep(NA, n_out)
      }
      pop <- rep(NA, n_out)
      
      count <- 0
      for(i in 1:n_subset){
        poly_i <- poly_use[i]
        which_times <- inc_times[inc_poly == poly_i]
        
        polygon_start_index[i] <- count
        n_match <- poly_size[i] * poly_times[i]
        which_i <- count + (1:n_match)
        count <- count + n_match
        polygon_end_index[i] <- count
        
        
        
        which_inc_match <- which(inc_poly == poly_use[i])
        times_match <- inc_times[which_inc_match]
        which_static_match <- which(poly_static == poly_use[i])
        which_dynamic_match <- which(poly_dynamic == poly_use[i])
        times_dynamic_match <- times_dynamic[which_dynamic_match]
        
        
        # if(i == 1) print(which_match)
        
        inc_vec <- rep(NA, n_match)
        for(j in 1:n_inc){
          for(t in 1:length(which_times)){
            t_i <- which_times[t]
            which_match <- which_inc_match[times_match == t_i]
            
            # which_match <- which(inc_poly == poly_use[i] & inc_times == t_i)
            inc_vec[(1:poly_size[i]) + (t-1)*poly_size[i]] <- rep(obs_list[[which_inc[j]]][which_match], poly_size[i])
          }
          out_list[[which_inc[j]]][which_i] <- inc_vec
        }
        
        pop_vec <- rep(NA, n_match)
        for(j in 1:n_dynamic){
          for(t in 1:length(which_times)){
            t_i <- which_times[t]
            # which_match <- which(poly_dynamic == poly_use[i] & times_dynamic == t_i)
            which_match <- which_dynamic_match[times_dynamic_match == t_i]
            
            inc_vec[(1:poly_size[i]) + (t-1)*poly_size[i]] <- obs_list[[which_dynamic[j]]][which_match]
            pop_vec[(1:poly_size[i]) + (t-1)*poly_size[i]] <- pop_dynamic[which_match]
          }
          out_list[[which_dynamic[j]]][which_i] <- inc_vec
        }
        
        
        for(j in 1:n_static){
          out_list[[which_static[j]]][which_i] <- rep(obs_list[[which_static[j]]][which_static_match],
                                                      each=poly_times[i])
        }
        
        pop[which_i] <- pop_vec
      }
    }
  }
  return(list(obs=out_list, pop=pop, index_use=index_use,
              start_index=polygon_start_index,
              end_index = polygon_end_index,
              pop=pop)
  )
}