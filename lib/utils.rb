require 'csv'
require 'tempfile'


module OpenTox

  module Algorithm

    @ambit_descriptor_algorithm_uri = "http://apps.ideaconsult.net:8080/ambit2/algorithm/org.openscience.cdk.qsar.descriptors.molecular."
    @ambit_ds_service_uri = "http://apps.ideaconsult.net:8080/ambit2/dataset/"
    @ambit_mopac_model_uri = "http://apps.ideaconsult.net:8080/ambit2/model/69632"
    @keysfile = File.join(ENV['HOME'], ".opentox", "config", "pc_descriptors.yaml")

    include OpenTox

    # Calculate physico-chemical descriptors.
    # @param[Hash] required: :dataset_uri, :pc_type, :rjb, :task, :add_uri, optional: :descriptor, :lib, :subjectid
    # @return[String] dataset uri
    def self.pc_descriptors(params)

      ds = OpenTox::Dataset.find(params[:dataset_uri],params[:subjectid])
      compounds = ds.compounds.collect
      task_weights = {"joelib"=> 20, "openbabel"=> 1, "cdk"=> 50 }
      task_weights.keys.each { |step| task_weights.delete(step) if (params[:lib] && (!params[:lib].split(",").include?(step)))}
      task_weights["load"] = 10
      task_sum = Float task_weights.values.sum
      task_weights.keys.each { |step| task_weights[step] /= task_sum }
      task_weights.keys.each { |step| task_weights[step] = (task_weights[step]*100).floor }
      
      jl_master=nil
      cdk_master=nil
      ob_master=nil


      # # # openbabel (via ruby bindings)
      if !params[:lib] || params[:lib].split(",").include?("openbabel")
        ob_master, ob_ids = get_ob_descriptors( { :compounds => compounds, :pc_type => params[:pc_type], :descriptor => params[:descriptor] } ) 
        params[:task].progress(params[:task].metadata[OT.percentageCompleted] + task_weights["openbabel"]) if params[:task]
      end


      # # # joelib (via rjb)
      if !params[:lib] || params[:lib].split(",").include?("joelib")
        jl_master, jl_ids = get_jl_descriptors( { :compounds => compounds, :rjb => params[:rjb][:jl], :pc_type => params[:pc_type], :descriptor => params[:descriptor] } ) 
        params[:task].progress(params[:task].metadata[OT.percentageCompleted] + task_weights["joelib"]) if params[:task]
      end


      # # # cdk (via rjb)
      if !params[:lib] || params[:lib].split(",").include?("cdk")
        cdk_master, cdk_ids = get_cdk_descriptors( { :compounds => compounds, :rjb => params[:rjb][:cdk], :pc_type => params[:pc_type], :descriptor => params[:descriptor] } ) 
        cdk_ids = cdk_master[0].collect { |id| id.to_s.split("-").first } # get column headers
        cdk_single_ids = cdk_master[0].collect { |id| id.to_s.sub(/[^-]*-/,"").gsub(/[\/.\\\(\)\{\}\[\]]/,"_") } # get column headers w/ nice '_'
        cdk_master[0] = cdk_single_ids.collect{|x| x} # Single IDs as features in result ds
        cdk_ids.shift # remove SMILES
        cdk_single_ids.shift # remove SMILES
        params[:task].progress(params[:task].metadata[OT.percentageCompleted] + task_weights["cdk"]) if params[:task]
      end

      # # # fuse CSVs ("master" structures)
      if jl_master && cdk_master
        nr_cols = (jl_master[0].size)-1
        LOGGER.debug "Merging #{nr_cols} new columns"
        cdk_master.each {|row| nr_cols.times { row.push(nil) }  }
        jl_master.each do |row|
          temp = cdk_master.assoc(row[0]) # Finds the appropriate line in master
          ((-1*nr_cols)..-1).collect.each { |idx|
            temp[idx] = row[nr_cols+idx+1] if temp # Updates columns if line is found
          }
        end
        master = cdk_master
      else # either jl_master or cdk_master nil
        master = jl_master || cdk_master
      end

      if ob_master && master
        nr_cols = (ob_master[0].size)-1
        LOGGER.debug "Merging #{nr_cols} new columns"
        master.each {|row| nr_cols.times { row.push(nil) }  } # Adds empty columns to all rows
        ob_master.each do |row|
          temp = master.assoc(row[0]) # Finds the appropriate line in master
          ((-1*nr_cols)..-1).collect.each { |idx|
            temp[idx] = row[nr_cols+idx+1] if temp # Updates columns if line is found
          }
        end
      else # either ob_master or master nil
        master = ob_master || master
      end

      if master

        ds = OpenTox::Dataset.find( 
          OpenTox::RestClientWrapper.post(
            File.join(CONFIG[:services]["opentox-dataset"]), master.collect { |row| row.join(",") }.join("\n"), {:content_type => "text/csv", :subjectid => params[:subjectid]}
          ),params[:subjectid]
        ) 

        # # # add feature metadata
        pc_descriptors = YAML::load_file(@keysfile)
        cdk_single_ids && cdk_single_ids.each_with_index { |id,idx|
          raise "Feature not found" if ! ds.features[File.join(ds.uri, "feature", id.to_s)]
          ds.add_feature_metadata(File.join(ds.uri, "feature", id.to_s),{DC.description => "#{pc_descriptors[cdk_ids[idx]][:name]} [#{pc_descriptors[cdk_ids[idx]][:pc_type]}, #{pc_descriptors[cdk_ids[idx]][:lib]}]"})
          creator_uri = ds.uri.gsub(/\/dataset\/.*/, "/algorithm/pc")
          creator_uri += "/#{cdk_ids[idx]}" if params[:add_uri]
          ds.add_feature_metadata(File.join(ds.uri, "feature", id.to_s),{DC.creator => creator_uri})
          ds.add_feature_metadata(File.join(ds.uri, "feature", id.to_s),{OT.hasSource => params[:dataset_uri]})
        }
        ob_ids && ob_ids.each { |id|
          raise "Feature not found" if ! ds.features[File.join(ds.uri, "feature", id.to_s)]
          ds.add_feature_metadata(File.join(ds.uri, "feature", id.to_s),{DC.description => "#{pc_descriptors[id][:name]} [#{pc_descriptors[id][:pc_type]}, #{pc_descriptors[id][:lib]}]"})
          creator_uri = ds.uri.gsub(/\/dataset\/.*/, "/algorithm/pc")
          creator_uri += "/#{id}" if params[:add_uri]
          ds.add_feature_metadata(File.join(ds.uri, "feature", id.to_s),{DC.creator => creator_uri})
          ds.add_feature_metadata(File.join(ds.uri, "feature", id.to_s),{OT.hasSource => params[:dataset_uri]})
        }
        jl_ids && jl_ids.each { |id|
          raise "Feature not found" if ! ds.features[File.join(ds.uri, "feature", id.to_s)]
          ds.add_feature_metadata(File.join(ds.uri, "feature", id.to_s),{DC.description => "#{pc_descriptors[id][:name]} [#{pc_descriptors[id][:pc_type]}, #{pc_descriptors[id][:lib]}]"})
          creator_uri = ds.uri.gsub(/\/dataset\/.*/, "/algorithm/pc")
          creator_uri += "/#{id}" if params[:add_uri]
          ds.add_feature_metadata(File.join(ds.uri, "feature", id.to_s),{DC.creator => creator_uri})
          ds.add_feature_metadata(File.join(ds.uri, "feature", id.to_s),{OT.hasSource => params[:dataset_uri]})
        }

        ds.save(params[:subjectid])
      else
        raise OpenTox::BadRequestError.new "No descriptors matching your criteria found."
      end

    end
    

    # Calculate OpenBabel physico-chemical descriptors.
    # @param[Hash] required: :compounds, :pc_type, :task, optional: :descriptor
    # @return[Array] CSV, array of field ids, array of field descriptions
    def self.get_ob_descriptors(params)

      master = nil

      begin
        csvfile = Tempfile.open(['ob_descriptors-','.csv'])

        pc_descriptors = YAML::load_file(@keysfile)
        ids = pc_descriptors.collect{ |id, info| 
          id if info[:lib] == "openbabel" && params[:pc_type].split(",").include?(info[:pc_type]) && (!params[:descriptor] || id == params[:descriptor])
        }.compact

        if ids.length > 0
          csvfile.puts((["SMILES"] + ids).join(","))
          
          # remember inchis
          inchis = params[:compounds].collect { |c_uri| 
            URI.encode_www_form_component(OpenTox::Compound.new(c_uri).to_inchi)
          }

          # Process compounds
          obmol = OpenBabel::OBMol.new
          obconversion = OpenBabel::OBConversion.new
          obconversion.set_in_and_out_formats 'inchi', 'can'

          inchis.each_with_index { |inchi, c_idx| 
            row = [inchis[c_idx]]
            obconversion.read_string(obmol, URI.decode_www_form_component(inchi))
            ids.each { |name|
              if obmol.respond_to?(name.underscore)
                val = eval("obmol.#{name.underscore}") if obmol.respond_to?(name.underscore) 
              else
                if name != "nF" && name != "spinMult" && name != "nHal" && name != "logP"
                  val = OpenBabel::OBDescriptor.find_type(name.underscore).predict(obmol)
                elsif name == "nF"
                  val = OpenBabel::OBDescriptor.find_type("nf").predict(obmol)
                elsif name == "spinMult" || name == "nHal" || name == "logP"
                  val = OpenBabel::OBDescriptor.find_type(name).predict(obmol)
                end
              end
              if OpenTox::Algorithm.numeric?(val)
                val = Float(val)
                val = nil if val.nan?
                val = nil if (val && val.infinite?)
              end
              row << val
            }
            LOGGER.debug "Compound #{c_idx+1} (#{inchis.size}), #{row.size} entries"
            csvfile.puts(row.join(","))
            csvfile.flush
          }
          master = CSV::parse(File.open(csvfile.path, "rb").read)
        end

      rescue Exception => e
        LOGGER.debug "#{e.class}: #{e.message}"
        LOGGER.debug "Backtrace:\n\t#{e.backtrace.join("\n\t")}"
      ensure
        csvfile.close!
      end

      [ master, ids ]

    end



    # Calculate Joelib2 physico-chemical descriptors.
    # @param[Hash] required: :compounds, :pc_type, :task, optional: :descriptor
    # @return[Array] CSV, array of field ids, array of field descriptions
    def self.get_jl_descriptors(params)

      master = nil
      s = params[:rjb]; raise "No Java environment" unless s

      # Load keys, enter CSV headers
      begin
        csvfile = Tempfile.open(['jl_descriptors-','.csv'])

        pc_descriptors = YAML::load_file(@keysfile)
        ids = pc_descriptors.collect{ |id, info| 
          id if info[:lib] == "joelib" && params[:pc_type].split(",").include?(info[:pc_type]) && (!params[:descriptor] || id == params[:descriptor])
        }.compact


        if ids.length > 0
          csvfile.puts((["SMILES"] + ids).join(","))

          # remember inchis
          inchis = params[:compounds].collect { |c_uri| 
            cmpd = OpenTox::Compound.new(c_uri)
            URI.encode_www_form_component(cmpd.to_inchi)
          }

          # Process compounds
          params[:compounds].each_with_index { |c_uri, c_idx| 
            cmpd = OpenTox::Compound.new(c_uri)
            inchi = cmpd.to_inchi
            sdf_data = cmpd.to_sdf

            infile = Tempfile.open(['jl_descriptors-in-','.sdf'])
            outfile_path = infile.path.gsub(/jl_descriptors-in/,"jl_descriptors-out")

            begin
              infile.puts sdf_data
              infile.flush
              s.new(infile.path, outfile_path) # runs joelib
                    
              row = [inchis[c_idx]]
              ids.each_with_index do |k,i| # Fill row
                re = Regexp.new(k)
                open(outfile_path) do |f|
                  f.each do |line|
                    if @prev == k
                      entry = line.chomp
                      val = nil
                      if OpenTox::Algorithm.numeric?(entry)
                        val = Float(entry)
                        val = nil if val.nan?
                        val = nil if (val && val.infinite?)
                      end
                      row << val
                      break
                    end
                    @prev = line.gsub(/^.*types./,"").gsub(/count./,"").gsub(/>/,"").chomp if line =~ re
                  end
                end
              end
              LOGGER.debug "Compound #{c_idx+1} (#{inchis.size}), #{row.size} entries"
              csvfile.puts(row.join(","))
              csvfile.flush

            rescue Exception => e
              LOGGER.debug "#{e.class}: #{e.message}"
              LOGGER.debug "Backtrace:\n\t#{e.backtrace.join("\n\t")}"
            ensure
              File.delete(infile.path.gsub(/\.sdf/,".numeric.sdf"))
              File.delete(outfile_path)
              infile.close!
            end
          }
          master = CSV::parse(File.open(csvfile.path, "rb").read)
        end

      rescue Exception => e
        LOGGER.debug "#{e.class}: #{e.message}"
        LOGGER.debug "Backtrace:\n\t#{e.backtrace.join("\n\t")}"
      ensure
        [ csvfile].each { |f| f.close! }
      end

      [ master, ids ]

    end

    
    def self.get_cdk_descriptors(params)

      master = nil
      t = params[:rjb]; raise "No Java environment" unless t

      # Load keys, enter CSV headers
      begin
        pc_descriptors = YAML::load_file(@keysfile)
        ids = pc_descriptors.collect{ |id, info| 
          id if info[:lib] == "cdk" && params[:pc_type].split(",").include?(info[:pc_type]) && (!params[:descriptor] || id == params[:descriptor])
        }.compact

        if ids.length > 0
          # remember inchis
          inchis = params[:compounds].collect { |c_uri| 
            cmpd = OpenTox::Compound.new(c_uri)
            URI.encode_www_form_component(cmpd.to_inchi)
          }
          begin
            # Process compounds
            sdf_data = []
            params[:compounds].each_with_index { |c_uri, c_idx| 
              cmpd = OpenTox::Compound.new(c_uri)
              LOGGER.debug "3D for #{cmpd.to_smiles}"
              obconv = OpenBabel::OBConversion.new
              obmol = OpenBabel::OBMol.new
              obconv.set_in_format("smi") 
              obconv.read_string(obmol, cmpd.to_smiles) 
              obconv.set_out_format("sdf") 
              sdf_string = obconv.write_string(obmol)  
              gen3d = OpenBabel::OBOp.find_type("Gen3D") 
              gen3d.do(obmol) 
              sdf_string_3d = obconv.write_string(obmol)  
              if sdf_string_3d.index(/.nan/).nil?
                sdf_data << sdf_string_3d
              else
                sdf_data << sdf_string
                LOGGER.debug "3D failed (using 2D)"
              end
            }
            infile = Tempfile.open(['jl_descriptors-in-','.sdf'])
            csvfile = infile.path.gsub(/jl_descriptors-in/,"jl_descriptors-out").gsub(/\.sdf/,".csv")
            infile.puts sdf_data.join("")
            infile.flush
            t.new(infile.path, csvfile, ids.join(',')) # runs cdk
          rescue Exception => e
            LOGGER.debug "#{e.class}: #{e.message}"
            LOGGER.debug "Backtrace:\n\t#{e.backtrace.join("\n\t")}"
          ensure
            infile.close!
          end
          master = CSV::parse(File.open(csvfile, "rb").read)
          master.each_with_index { |row, idx|
            if idx != 0 # not alter headers
              row[0] = inchis[idx-1]
              row.collect! { |x| x.to_s == "null" ? nil : x }
            end
          }
        end
      rescue Exception => e
        LOGGER.debug "#{e.class}: #{e.message}"
        LOGGER.debug "Backtrace:\n\t#{e.backtrace.join("\n\t")}"
      ensure
        [ csvfile ].each { |f| File.delete(f) }
      end

      [ master, ids ]

    end


    # Gauss kernel
    # @return [Float] 
    def self.gauss(x, sigma = 0.3) 
      d = 1.0 - x.to_f
      Math.exp(-(d*d)/(2*sigma*sigma))
    end


    # For symbolic features
    # @param [Array] Array to test, must indicate non-occurrence with 0.
    # @return [Boolean] Whether the feature is singular or non-occurring or present everywhere.
    def self.isnull_or_singular?(array)
      nr_zeroes = array.count(0)
      return (nr_zeroes == array.size) ||    # remove non-occurring feature
             (nr_zeroes == array.size-1) ||  # remove singular feature
             (nr_zeroes == 0)                # also remove feature present everywhere
    end


    # Numeric value test
    # @param[Object] value
    # @return [Boolean] Whether value is a number
    def self.numeric?(value)
      true if Float(value) rescue false
    end


    # For symbolic features
    # @param [Array] Array to test, must indicate non-occurrence with 0.
    # @return [Boolean] Whether the feature has variance zero.
    def self.zero_variance?(array)
      return array.uniq.size == 1
    end
    

    # Sum of an array for Arrays.
    # @param [Array] Array with values
    # @return [Integer] Sum of size of values
    def self.sum_size(array)
      sum=0
      array.each { |e| sum += e.size }
      return sum
    end


    # Minimum Frequency
    # @param [Integer] per-mil value
    # return [Integer] min-frequency
    def self.min_frequency(training_dataset,per_mil)
      minfreq = per_mil * training_dataset.compounds.size.to_f / 1000.0 # AM sugg. 8-10 per mil for BBRC, 50 per mil for LAST
      minfreq = 2 unless minfreq > 2
      Integer (minfreq)
    end


    # Effect calculation for classification. It is assumed that the elements of the arrays match each other pairwise
    # @param [Array] Array of occurrences per class (in the form of Enumerables).
    # @param [Array] Array of database instance counts per class.
    def self.effect(occurrences, db_instances)
      max=0
      max_value=0
      nr_o = self.sum_size(occurrences)
      nr_db = db_instances.to_scale.sum

      occurrences.each_with_index { |o,i| # fminer outputs occurrences sorted reverse by activity.
        actual = o.size.to_f/nr_o
        expected = db_instances[i].to_f/nr_db
        if actual > expected
          if ((actual - expected) / actual) > max_value
           max_value = (actual - expected) / actual # 'Schleppzeiger'
            max = i
          end
        end
      }
      max
    end
    

    # neighbors

    module Neighbors
      
      # Get confidence.
      # @param[Hash] Required keys: :sims, :acts
      # @return[Float] Confidence
      def self.get_confidence(params)
        conf = params[:sims].inject{|sum,x| sum + x }
        confidence = conf/params[:sims].size
        LOGGER.debug "Confidence is: '" + confidence.to_s + "'."
        return confidence
      end

    end


    # Similarity calculations
    module Similarity

      # Tanimoto similarity
      # @param [Hash, Array] fingerprints of first compound
      # @param [Hash, Array] fingerprints of second compound
      # @return [Float] (Weighted) tanimoto similarity
      def self.tanimoto(fingerprints_a,fingerprints_b,weights=nil,params=nil)

        common_p_sum = 0.0
        all_p_sum = 0.0

        # fingerprints are hashes
        if fingerprints_a.class == Hash && fingerprints_b.class == Hash
          common_features = fingerprints_a.keys & fingerprints_b.keys
          all_features = (fingerprints_a.keys + fingerprints_b.keys).uniq
          if common_features.size > 0
            common_features.each{ |f| common_p_sum += [ fingerprints_a[f], fingerprints_b[f] ].min }
            all_features.each{ |f| all_p_sum += [ fingerprints_a[f],fingerprints_b[f] ].compact.max } # compact, since one fp may be empty at that pos
          end

        # fingerprints are arrays
        elsif fingerprints_a.class == Array && fingerprints_b.class == Array
          size = [ fingerprints_a.size, fingerprints_b.size ].min
          LOGGER.warn "fingerprints don't have equal size" if fingerprints_a.size != fingerprints_b.size
          (0...size).each { |idx|
            common_p_sum += [ fingerprints_a[idx], fingerprints_b[idx] ].min
            all_p_sum += [ fingerprints_a[idx], fingerprints_b[idx] ].max
          }
        end

        (all_p_sum > 0.0) ? (common_p_sum/all_p_sum) : 0.0

      end


      # Cosine similarity
      # @param [Hash] properties_a key-value properties of first compound
      # @param [Hash] properties_b key-value properties of second compound
      # @return [Float] cosine of angle enclosed between vectors induced by keys present in both a and b
      def self.cosine(fingerprints_a,fingerprints_b,weights=nil)

        # fingerprints are hashes
        if fingerprints_a.class == Hash && fingerprints_b.class == Hash
          a = []; b = []
          common_features = fingerprints_a.keys & fingerprints_b.keys
          if common_features.size > 1
            common_features.each do |p|
              a << fingerprints_a[p]
              b << fingerprints_b[p]
            end
          end

        # fingerprints are arrays
        elsif fingerprints_a.class == Array && fingerprints_b.class == Array
          a = fingerprints_a
          b = fingerprints_b
        end

        (a.size > 0 && b.size > 0) ? self.cosine_num(a.to_gv, b.to_gv) : 0.0

      end


      # Cosine similarity
      # @param [GSL::Vector] a
      # @param [GSL::Vector] b
      # @return [Float] cosine of angle enclosed between a and b
      def self.cosine_num(a, b)
        if a.size>12 && b.size>12
          a = a[0..11]
          b = b[0..11]
        end
        a.dot(b) / (a.norm * b.norm)
      end


      # Outlier detection based on Mahalanobis distances
      # Multivariate detection on X, univariate detection on y
      # Uses an existing Rinruby instance, if possible
      # @param[Hash] Keys query_matrix, data_matrix, acts are required; r, p_outlier optional
      # @return[Array] indices identifying outliers (may occur several times, this is intended)
      def self.outliers(params)
        outlier_array = []
        data_matrix = params[:data_matrix]
        query_matrix = params[:query_matrix]
        acts = params[:acts]
        begin
          LOGGER.debug "Outliers (p=#{params[:p_outlier] || 0.9999})..."
          r = ( params[:r] || RinRuby.new(false,false) )
          r.eval "suppressPackageStartupMessages(library(\"robustbase\"))"
          r.eval "outlier_threshold = #{params[:p_outlier] || 0.999}"
          nr_cases, nr_features = data_matrix.to_a.size, data_matrix.to_a[0].size
          r.odx = data_matrix.to_a.flatten
          r.q = query_matrix.to_a.flatten
          r.y = acts.to_a.flatten
          r.eval "odx = matrix(odx, #{nr_cases}, #{nr_features}, byrow=T)"
          r.eval 'odx = rbind(q,odx)' # query is nr 0 (1) in ruby (R)
          r.eval 'mah = covMcd(odx)$mah' # run MCD alg
          r.eval "mah = pchisq(mah,#{nr_features})"
          r.eval 'outlier_array = which(mah>outlier_threshold)'  # multivariate outliers using robust mahalanobis
          outlier_array = r.outlier_array.to_a.collect{|v| v-2 }  # translate to ruby index (-1 for q, -1 due to ruby)
          r.eval 'fqu = matrix(summary(y))[2]'
          r.eval 'tqu = matrix(summary(y))[5]'
          r.eval 'outlier_array = which(y>(tqu+1.5*IQR(y)))'     # univariate outliers due to Tukey (http://goo.gl/mwzNH)
          outlier_array += r.outlier_array.to_a.collect{|v| v-1 } # translate to ruby index (-1 due to ruby)
          r.eval 'outlier_array = which(y<(fqu-1.5*IQR(y)))'
          outlier_array += r.outlier_array.to_a.collect{|v| v-1 }
        rescue Exception => e
          LOGGER.debug "#{e.class}: #{e.message}"
          #LOGGER.debug "Backtrace:\n\t#{e.backtrace.join("\n\t")}"
        end
        outlier_array
      end

    end


  end

end

