module OpenTox
    module Transform
    
      # LogAutoScaler for GSL vectors.
      # Take log and scale.
      # Uses Statsample Library (http://ruby-statsample.rubyforge.org/) by C. Bustos
      class LogAutoScale
        attr_accessor :vs, :offset, :autoscaler

        # @param [GSL::Vector] Values to transform using LogAutoScaling.
        def initialize values
          @distance_to_zero = 1.0
          begin
            raise "Cannot transform, values empty." if values.size==0
            vs = values.clone
            @offset = vs.min - @distance_to_zero
            @autoscaler = OpenTox::Transform::AutoScale.new mvlog(vs)
            @vs = @autoscaler.vs
          rescue Exception => e
            LOGGER.debug "#{e.class}: #{e.message}"
            LOGGER.debug "Backtrace:\n\t#{e.backtrace.join("\n\t")}"
          end
        end

        # @param [GSL::Vector] values to restore.
        # @return [GSL::Vector] transformed values.
        def restore values
          begin
            raise "Cannot transform, values empty." if values.size==0
            vs = values.clone
            rv = @autoscaler.restore(vs)
            rv.to_a.collect { |v| (10**v) + @offset }.to_gv
          rescue Exception => e
            LOGGER.debug "#{e.class}: #{e.message}"
            LOGGER.debug "Backtrace:\n\t#{e.backtrace.join("\n\t")}"
          end
        end

        # @param [GSL::Vector] values to transform.
        # @return [GSL::Vector] transformed values.
        def mvlog values 
          values.to_a.collect { |v| Math::log10(v - @offset) }.to_gv
        end

      end


      # Auto-Scaler for GSL vectors.
      # Center on mean and divide by standard deviation.
      # Uses Statsample Library (http://ruby-statsample.rubyforge.org/) by C. Bustos
      class AutoScale 
        attr_accessor :vs, :mean, :stdev

        # @param [GSL::Vector] values to transform using AutoScaling.
        def initialize values
          begin
            raise "Cannot transform, values empty." if values.size==0
            vs = values.clone
            @mean = vs.to_scale.mean
            @stdev = vs.to_scale.standard_deviation_population
            @stdev = 0.0 if @stdev.nan? 
            @vs = transform vs
          rescue Exception => e
            LOGGER.debug "#{e.class}: #{e.message}"
            LOGGER.debug "Backtrace:\n\t#{e.backtrace.join("\n\t")}"
          end
        end

        # @param [GSL::Vector] values to transform.
        # @return [GSL::Vector] transformed values.
        def transform values
          begin
            raise "Cannot transform, values empty." if values.size==0
            autoscale values.clone
          rescue Exception => e
            LOGGER.debug "#{e.class}: #{e.message}"
            LOGGER.debug "Backtrace:\n\t#{e.backtrace.join("\n\t")}"
          end
        end

        # @param [GSL::Vector] Values to restore.
        # @return [GSL::Vector] transformed values.
        def restore values
          begin
            raise "Cannot transform, values empty." if values.size==0
            rv_ss = values.clone.to_scale * @stdev unless @stdev == 0.0
            (rv_ss + @mean).to_gsl
          rescue Exception => e
            LOGGER.debug "#{e.class}: #{e.message}"
            LOGGER.debug "Backtrace:\n\t#{e.backtrace.join("\n\t")}"
          end
        end

        # @param [GSL::Vector] values to transform.
        # @return [GSL::Vector] transformed values.
        def autoscale values
          vs_ss = values.clone.to_scale - @mean
          @stdev == 0.0 ? vs_ss.to_gsl : ( vs_ss * ( 1 / @stdev) ).to_gsl
        end

      end


      # Principal Components Analysis.
      # Uses Statsample Library (http://ruby-statsample.rubyforge.org/) by C. Bustos
      class PCA
        attr_accessor :data_matrix, :data_transformed_matrix, :eigenvector_matrix, :eigenvalue_sums, :autoscaler

        # Creates a transformed dataset as GSL::Matrix.
        #
        # @param [GSL::Matrix] Data matrix.
        # @param [Float] Compression ratio from [0,1], default 0.05.
        # @return [GSL::Matrix] Data transformed matrix.
        def initialize data_matrix, compression=0.05, maxcols=(1.0/0.0)
          begin
            @data_matrix = data_matrix.clone
            @compression = compression.to_f
            @mean = Array.new
            @autoscaler = Array.new
            @cols = Array.new
            @maxcols = maxcols

            # Objective Feature Selection
            raise "Error! PCA needs at least two dimensions." if data_matrix.size2 < 2
            @data_matrix_selected = nil
            (0..@data_matrix.size2-1).each { |i|
              if !Algorithm::zero_variance?(@data_matrix.col(i).to_a)
                if @data_matrix_selected.nil?
                  @data_matrix_selected = GSL::Matrix.alloc(@data_matrix.size1, 1) 
                  @data_matrix_selected.col(0)[0..@data_matrix.size1-1] = @data_matrix.col(i)
                else
                  @data_matrix_selected = @data_matrix_selected.horzcat(GSL::Matrix.alloc(@data_matrix.col(i).to_a,@data_matrix.size1, 1))
                end
                @cols << i
              end             
            }
            raise "Error! PCA needs at least two dimensions." if (@data_matrix_selected.nil? || @data_matrix_selected.size2 < 2)

            # PCA uses internal centering on 0
            @data_matrix_scaled = GSL::Matrix.alloc(@data_matrix_selected.size1, @cols.size)
            (0..@cols.size-1).each { |i|
              as = OpenTox::Transform::AutoScale.new(@data_matrix_selected.col(i))
              @data_matrix_scaled.col(i)[0..@data_matrix.size1-1] = as.vs * as.stdev # re-adjust by stdev
              @mean << as.mean
              @autoscaler << as
            }

            # PCA
            data_matrix_hash = Hash.new
            (0..@cols.size-1).each { |i|
              column_view = @data_matrix_scaled.col(i)
              data_matrix_hash[i] = column_view.to_scale
            }
            dataset_hash = data_matrix_hash.to_dataset # see http://goo.gl/7XcW9
            cor_matrix=Statsample::Bivariate.correlation_matrix(dataset_hash)
            pca=Statsample::Factor::PCA.new(cor_matrix)

            # Select best eigenvectors
            pca.eigenvalues.each { |ev| raise "PCA failed!" unless !ev.nan? }
            @eigenvalue_sums = Array.new
            (0..@cols.size-1).each { |i|
              @eigenvalue_sums << pca.eigenvalues[0..i].inject{ |sum, ev| sum + ev }
            }
            eigenvectors_selected = Array.new
            pca.eigenvectors.each_with_index { |ev, i|
              if (@eigenvalue_sums[i] <= ((1.0-@compression)*@cols.size)) || (eigenvectors_selected.size == 0)
                eigenvectors_selected << ev.to_a unless @maxcols <= eigenvectors_selected.size
              end
            }
            @eigenvector_matrix = GSL::Matrix.alloc(eigenvectors_selected.flatten, eigenvectors_selected.size, @cols.size).transpose
            @data_transformed_matrix = (@eigenvector_matrix.transpose * @data_matrix_scaled.transpose).transpose

          rescue Exception => e
              LOGGER.debug "#{e.class}: #{e.message}"
              LOGGER.debug "Backtrace:\n\t#{e.backtrace.join("\n\t")}"
          end
        end

        # Transforms data to feature space found by PCA.
        #
        # @param [GSL::Matrix] Data matrix.
        # @return [GSL::Matrix] Transformed data matrix.
        def transform values
          begin
            vs = values.clone
            raise "Error! Too few columns for transformation." if vs.size2 < @cols.max
            data_matrix_scaled = GSL::Matrix.alloc(vs.size1, @cols.size)
            @cols.each_with_index { |i,j|
              data_matrix_scaled.col(j)[0..data_matrix_scaled.size1-1] = @autoscaler[j].transform(vs.col(i).to_a) * @autoscaler[j].stdev
            }
            (@eigenvector_matrix.transpose * data_matrix_scaled.transpose).transpose
          rescue Exception => e
            LOGGER.debug "#{e.class}: #{e.message}"
            LOGGER.debug "Backtrace:\n\t#{e.backtrace.join("\n\t")}"
          end
        end

        # Restores data in the original feature space (possibly with compression loss).
        #
        # @param [GSL::Matrix] Transformed data matrix.
        # @return [GSL::Matrix] Data matrix.
        def restore
          begin 
            data_matrix_restored = (@eigenvector_matrix * @data_transformed_matrix.transpose).transpose # reverse pca
            # reverse scaling
            (0..@cols.size-1).each { |i|
              data_matrix_restored.col(i)[0..data_matrix_restored.size1-1] += @mean[i]
            }
            data_matrix_restored
          rescue Exception => e
            LOGGER.debug "#{e.class}: #{e.message}"
            LOGGER.debug "Backtrace:\n\t#{e.backtrace.join("\n\t")}"
          end
        end

      end


      # Singular Value Decomposition
      # Uses Statsample Library (http://ruby-statsample.rubyforge.org/) by C. Bustos
      class SVD
        attr_accessor :data_matrix, :compression, :data_transformed_matrix, :uk, :vk, :eigk, :eigk_inv

        # Creates a transformed dataset as GSL::Matrix.
        #
        # @param [GSL::Matrix] Data matrix
        # @param [Float] Compression ratio from [0,1], default 0.20
        # @return [GSL::Matrix] Data transformed matrix

        def initialize data_matrix, compression=0.20
          begin
            @data_matrix = data_matrix.clone
            @compression = compression

            # Compute the SV Decomposition X=USV
            # vt is *not* the transpose of V here, but V itself (see http://goo.gl/mm2xz)!
            u, vt, s = data_matrix.SV_decomp 
            
            # Determine cutoff index
            s2 = s.mul(s) ; s2_sum = s2.sum
            s2_run = 0
            k = s2.size - 1
            s2.to_a.reverse.each { |v| 
              s2_run += v
              frac = s2_run / s2_sum
              break if frac > compression
              k -= 1
            }
            k += 1 if k == 0 # avoid uni-dimensional (always cos sim of 1)
            
            # Take the k-rank approximation of the Matrix
            #   - Take first k columns of u
            #   - Take first k columns of vt
            #   - Take the first k eigenvalues
            @uk = u.submatrix(nil, (0..k)) # used to transform column format data
            @vk = vt.submatrix(nil, (0..k)) # used to transform row format data
            s = GSL::Matrix.diagonal(s)
            @eigk = s.submatrix((0..k), (0..k))
            @eigk_inv = @eigk.inv

            # Transform data
            @data_transformed_matrix = @data_matrix * @vk * @eigk_inv # = @uk for all SVs

          rescue Exception => e
            LOGGER.debug "#{e.class}: #{e.message}"
            LOGGER.debug "Backtrace:\n\t#{e.backtrace.join("\n\t")}"
          end
        end


        # Transforms data instance (1 row) to feature space found by SVD.
        #
        # @param [GSL::Matrix] Data matrix (1 x m).
        # @return [GSL::Matrix] Transformed data matrix.
        def transform_instance values
          begin
            values * @vk * @eigk_inv
          rescue Exception => e
            LOGGER.debug "#{e.class}: #{e.message}"
            LOGGER.debug "Backtrace:\n\t#{e.backtrace.join("\n\t")}"
          end
        end
        alias :transform :transform_instance # make this the default (see PCA interface)

        # Transforms data feature (1 column) to feature space found by SVD.
        #
        # @param [GSL::Matrix] Data matrix (1 x n).
        # @return [GSL::Matrix] Transformed data matrix.
        def transform_feature values
          begin
            values * @uk * @eigk_inv
          rescue Exception => e
            LOGGER.debug "#{e.class}: #{e.message}"
            LOGGER.debug "Backtrace:\n\t#{e.backtrace.join("\n\t")}"
          end
        end


        # Restores data in the original feature space (possibly with compression loss).
        #
        # @param [GSL::Matrix] Transformed data matrix.
        # @return [GSL::Matrix] Data matrix.
        def restore
          begin 
            @data_transformed_matrix * @eigk * @vk.transpose  # reverse svd
          rescue Exception => e
            LOGGER.debug "#{e.class}: #{e.message}"
            LOGGER.debug "Backtrace:\n\t#{e.backtrace.join("\n\t")}"
          end
        end


      end

    end
end
