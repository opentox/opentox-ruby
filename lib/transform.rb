module OpenTox
    module Transform
    
      # LogAutoScaler for ruby Arrays.
      # Take log and scale.
      # Uses Statsample Library (http://ruby-statsample.rubyforge.org/) by C. Bustos
      class LogAutoScale
        attr_accessor :sv, :rv, :offset, :autoscaler

        # @params[Array] Values to transform using LogAutoScaling.
        def initialize values
          @distance_to_zero = 1.0
          begin
            raise "Cannot transform, values empty." if values.size==0
            @sv = values.to_scale
            @offset = @sv.min - @distance_to_zero
            @sv = @sv.to_a.collect { |v| Math::log10(v - @offset) }
            @autoscaler = OpenTox::Algorithm::Transform::AutoScale.new(@sv)
            @sv = @autoscaler.sv
          rescue Exception => e
            LOGGER.debug "#{e.class}: #{e.message}"
            LOGGER.debug "Backtrace:\n\t#{e.backtrace.join("\n\t")}"
          end
        end

        # @params[Array] ruby array of values to transform.
        # @returns[Array] ruby array of transformed values.
        def transform values
          values.collect! { |v| Math::log10(v - @offset) }
          @autoscaler.transform values
        end

        # @params[Array] ruby array of values to restore.
        # @returns[Array] ruby array of transformed values.
        def restore values
          @rv = @autoscaler.restore values
          @rv.collect! { |v| (10**v) + @offset }
        end

      end


      # Auto-Scaler for ruby Arrays.
      # Center on mean and divide by standard deviation.
      # Uses Statsample Library (http://ruby-statsample.rubyforge.org/) by C. Bustos
      class AutoScale 
        attr_accessor :sv, :rv, :mean, :stdev

        # @params[Array] ruby array of values to transform using AutoScaling.
        def initialize values
          begin
            raise "Cannot transform, values empty." if values.size==0
            @sv = values.to_scale
            @mean = @sv.mean
            @stdev = @sv.standard_deviation_population
            @sv = transform @sv.to_a
          rescue Exception => e
            LOGGER.debug "#{e.class}: #{e.message}"
            LOGGER.debug "Backtrace:\n\t#{e.backtrace.join("\n\t")}"
          end
        end

        # @params[Array] ruby array of values to transform.
        # @returns[Array] ruby array of transformed values.
        def transform values
          values = values.to_scale - @mean
          values = values * ( 1 / @stdev) unless @stdev == 0.0
          values.to_a
        end

        # @params[Array] Values to restore.
        # @returns[Array] ruby array of transformed values.
        def restore values
          @rv = values.to_scale
          @rv = @rv * @stdev unless @stdev == 0.0
          @rv = (@rv + @mean).to_a
        end
      end


      # Principal Components Analysis
      # Uses Statsample Library (http://ruby-statsample.rubyforge.org/) by C. Bustos
      class PCA
        attr_accessor :data_matrix, :data_transformed_matrix, :eigenvector_matrix, :eigenvalue_sums, :autoscaler

        # Creates a transformed dataset as GSL::Matrix.
        # @param [GSL::Matrix] Data matrix.
        # @param [Float] Compression ratio from [0,1].
        # @return [GSL::Matrix] Data transformed matrix.
        def initialize data_matrix, compression=0.05
          begin
            @data_matrix = data_matrix
            @compression = compression.to_f
            @stdev = Array.new
            @mean = Array.new

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
              end             
            }
            raise "Error! PCA needs at least two dimensions." if (@data_matrix_selected.nil? || @data_matrix_selected.size2 < 2)

            # Scaling of Axes
            @data_matrix_scaled = GSL::Matrix.alloc(@data_matrix_selected.size1, @data_matrix_selected.size2)
            (0..@data_matrix_selected.size2-1).each { |i|
              @autoscaler = OpenTox::Transform::AutoScale.new(@data_matrix_selected.col(i))
              @data_matrix_scaled.col(i)[0..@data_matrix.size1-1] = @autoscaler.sv
              @stdev << @autoscaler.stdev
              @mean << @autoscaler.mean
            }
            data_matrix_hash = Hash.new
            (0..@data_matrix_scaled.size2-1).each { |i|
              column_view = @data_matrix_scaled.col(i)
              data_matrix_hash[i] = column_view.to_scale
            }
            dataset_hash = data_matrix_hash.to_dataset # see http://goo.gl/7XcW9
            cor_matrix=Statsample::Bivariate.correlation_matrix(dataset_hash)
            pca=Statsample::Factor::PCA.new(cor_matrix)
            pca.eigenvalues.each { |ev| raise "PCA failed!" unless !ev.nan? }
            @eigenvalue_sums = Array.new
            (0..dataset_hash.fields.size-1).each { |i|
              @eigenvalue_sums << pca.eigenvalues[0..i].inject{ |sum, ev| sum + ev }
            }
            eigenvectors_selected = Array.new
            pca.eigenvectors.each_with_index { |ev, i|
              if (@eigenvalue_sums[i] <= ((1.0-@compression)*dataset_hash.fields.size)) || (eigenvectors_selected.size == 0)
                eigenvectors_selected << ev.to_a
              end
            }
            @eigenvector_matrix = GSL::Matrix.alloc(eigenvectors_selected.flatten, eigenvectors_selected.size, dataset_hash.fields.size).transpose
            dataset_matrix = dataset_hash.to_gsl.transpose
            @data_transformed_matrix = (@eigenvector_matrix.transpose * dataset_matrix).transpose

          rescue Exception => e
              LOGGER.debug "#{e.class}: #{e.message}"
              LOGGER.debug "Backtrace:\n\t#{e.backtrace.join("\n\t")}"
          end
        end

        def transform values
          data_matrix_scaled = GSL::Matrix.alloc(values.size1, values.size2)
          (0..values.size2-1).each { |i|
            autoscaler = OpenTox::Transform::AutoScale.new(values.col(i))
            data_matrix_scaled.col(i)[0..data_matrix_scaled.size1-1] = autoscaler.sv
          }
          (@eigenvector_matrix.transpose * data_matrix_scaled.transpose).transpose
        end

        # Restores data in the original feature space (possibly with compression loss).
        # @return [GSL::Matrix] Data matrix.
        def restore
          begin 
            data_matrix_restored = (@eigenvector_matrix * @data_transformed_matrix.transpose).transpose # reverse pca
            # reverse scaling
            (0..data_matrix_restored.size2-1).each { |i|
              data_matrix_restored.col(i)[0..data_matrix_restored.size1-1] *= @stdev[i] unless @stdev[i] == 0.0
              data_matrix_restored.col(i)[0..data_matrix_restored.size1-1] += @mean[i]
            }
            data_matrix_restored
          rescue Exception => e
            LOGGER.debug "#{e.class}: #{e.message}"
            LOGGER.debug "Backtrace:\n\t#{e.backtrace.join("\n\t")}"
          end
        end

      end

    end
end
