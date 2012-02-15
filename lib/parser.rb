require 'spreadsheet'
require 'roo'

class String

  # Split RDF statement into triples
  # @return [Array] Array with [subject,predicate,object]
  def to_triple
    self.chomp.split(' ',3).collect{|i| i.sub(/\s+.$/,'').gsub(/[<>"]/,'')}
  end

end

module OpenTox

  # Parser for various input formats
  module Parser

    # OWL-DL parser 
    module Owl

      # Create a new OWL-DL parser
      # @param uri URI of OpenTox object
      # @return [OpenTox::Parser::Owl] OWL-DL parser
      def initialize(uri)
        @uri = uri
        @metadata = {}
      end

      # Read metadata from opentox service
      # @return [Hash] Object metadata
      def load_metadata(subjectid=nil)
        # avoid using rapper directly because of 2 reasons:
        # * http errors wont be noticed
        # * subjectid cannot be sent as header
        ##uri += "?subjectid=#{CGI.escape(subjectid)}" if subjectid 
        ## `rapper -i rdfxml -o ntriples #{uri} 2>/dev/null`.each_line do |line|
        if File.exist?(@uri)
          file = File.new(@uri)
        else
          file = Tempfile.new("ot-rdfxml")
          if @dataset
            uri = URI::parse(@uri)
            #remove params like dataset/<id>?max=3 from uri, not needed for metadata
            uri.query = nil 
            uri.path = File.join(uri.path,"metadata")
            uri = uri.to_s
          else
            uri = @uri
          end
          file.puts OpenTox::RestClientWrapper.get uri,{:subjectid => subjectid,:accept => "application/rdf+xml"},nil,false
          file.close
          to_delete = file.path
        end
        statements = []
        parameter_ids = []
        `rapper -i rdfxml -o ntriples #{file.path} 2>/dev/null`.each_line do |line|
          triple = line.to_triple
          if triple[0] == @uri
            if triple[1] == RDF.type || triple[1]==OT.predictedVariables || triple[1]==OT.independentVariables # allow multiple types
              @metadata[triple[1]] = [] unless @metadata[triple[1]]
              @metadata[triple[1]] << triple[2].split('^^').first
            else
              @metadata[triple[1]] = triple[2].split('^^').first
            end
          end
          statements << triple 
          parameter_ids << triple[2] if triple[1] == OT.parameters
        end
        File.delete(to_delete) if to_delete
        unless parameter_ids.empty?
          @metadata[OT.parameters] = []
          parameter_ids.each do |p|
            parameter = {}
            statements.each{ |t| parameter[t[1]] = t[2] if t[0] == p and t[1] != RDF['type']}
            @metadata[OT.parameters] << parameter
          end
        end
        #@metadata.each do |k,v|
          #v = v.first if v and v.size == 1
        #end
        @metadata
      end
      
      # creates owl object from rdf-data
      # @param [String] rdf
      # @param [String] type of the info (e.g. OT.Task, OT.ErrorReport) needed to get the subject-uri
      # @return [Owl] with uri and metadata set 
      def self.from_rdf( rdf, type, allow_multiple = false )

        uris = Array.new
        owls = Array.new

        # write to file and read convert with rapper into tripples
        file = Tempfile.new("ot-rdfxml")
        file.puts rdf
        file.close
        #puts "cmd: rapper -i rdfxml -o ntriples #{file} 2>/dev/null"
        triples = `rapper -i rdfxml -o ntriples #{file.path} 2>/dev/null`
        
        # load uri via type
        uri = nil
        triples.each_line do |line|
          triple = line.to_triple
          if triple[1] == RDF['type'] and triple[2]==type
             if !allow_multiple
               raise "uri already set, two uris found with type: "+type.to_s if uri
             end
             uri = triple[0]
             uris << uri
          end
        end
        File.delete(file.path)

        # load metadata
        uris.each { |uri|
          metadata = {}
          triples.each_line do |line|
            triple = line.to_triple
            metadata[triple[1]] = triple[2].split('^^').first if triple[0] == uri and triple[1] != RDF['type']
          end
          owl = Owl::Generic.new(uri)
          owl.metadata = metadata
          owls << owl
        }
        allow_multiple ? owls : owls[0]
      end
      
      # Generic parser for all OpenTox classes
      class Generic
        include Owl
        
        attr_accessor :uri, :metadata
      end

      # OWL-DL parser for datasets
      class Dataset

        include Owl

        attr_writer :uri

        # Create a new OWL-DL dataset parser
        # @param uri Dataset URI 
        # @return [OpenTox::Parser::Owl::Dataset] OWL-DL parser
        def initialize(uri, subjectid=nil)
          super uri
          @dataset = ::OpenTox::Dataset.new(@uri, subjectid)
        end

        # Read data from dataset service. Files can be parsed by setting #uri to a filename (after initialization with a real URI)
        # @example Read data from an external service
        #   parser = OpenTox::Parser::Owl::Dataaset.new "http://wwbservices.in-silico.ch/dataset/1"
        #   dataset = parser.load_uri
        # @example Create dataset from RDF/XML file
        #   dataset = OpenTox::Dataset.create
        #   parser = OpenTox::Parser::Owl::Dataaset.new dataset.uri
        #   parser.uri = "dataset.rdfxml" # insert your input file
        #   dataset = parser.load_uri
        #   dataset.save
        # @return [Hash] Internal dataset representation
        def load_uri(subjectid=nil)
          
          # avoid using rapper directly because of 2 reasons:
          # * http errors wont be noticed
          # * subjectid cannot be sent as header
          ##uri += "?subjectid=#{CGI.escape(subjectid)}" if subjectid
          ##`rapper -i rdfxml -o ntriples #{file} 2>/dev/null`.each_line do |line| 
          if File.exist?(@uri)
            file = File.new(@uri)
          else
            file = Tempfile.new("ot-rdfxml")
            file.puts OpenTox::RestClientWrapper.get @uri,{:subjectid => subjectid,:accept => "application/rdf+xml"},nil,false
            file.close
            to_delete = file.path
          end
          
          data = {}
          feature_values = {}
          feature = {}
          feature_accept_values = {}
          other_statements = {}
          `rapper -i rdfxml -o ntriples #{file.path} 2>/dev/null`.each_line do |line|
            triple = line.chomp.split(' ',3)
            triple = triple[0..2].collect{|i| i.sub(/\s+.$/,'').gsub(/[<>"]/,'')}
            case triple[1] 
            when /#{OT.values}/i
              data[triple[0]] = {:compound => "", :values => []} unless data[triple[0]]
              data[triple[0]][:values] << triple[2]  
            when /#{OT.value}/i
              feature_values[triple[0]] = triple[2] 
            when /#{OT.compound}/i
              data[triple[0]] = {:compound => "", :values => []} unless data[triple[0]]
              data[triple[0]][:compound] = triple[2]  
            when /#{OT.feature}/i
              feature[triple[0]] = triple[2]
            when /#{RDF.type}/i
              if triple[2]=~/#{OT.Compound}/i and !data[triple[0]]
                data[triple[0]] = {:compound => triple[0], :values => []} 
              end
            when /#{OT.acceptValue}/i # acceptValue in ambit datasets is only provided in dataset/<id> no in dataset/<id>/features  
              feature_accept_values[triple[0]] = [] unless feature_accept_values[triple[0]]
              feature_accept_values[triple[0]] << triple[2]
            else 
            end
          end
          File.delete(to_delete) if to_delete
          data.each do |id,entry|
            if entry[:values].size==0
              # no feature values add plain compounds
              @dataset.add_compound(entry[:compound])
            else
              entry[:values].each do |value_id|
                if feature_values[value_id]
                  split = feature_values[value_id].split(/\^\^/)
                  case split[-1]
                  when XSD.double, XSD.float 
                    value = split.first.to_f
                  when XSD.boolean
                    value = split.first=~/(?i)true/ ? true : false                
                  else
                    value = split.first
                  end
                end
                @dataset.add entry[:compound],feature[value_id],value
              end
            end
          end
          load_features subjectid
          feature_accept_values.each do |feature, values|
            @dataset.features[feature][OT.acceptValue] = values
          end
          @dataset.metadata = load_metadata(subjectid)
          @dataset
        end

        # Read only features from a dataset service. 
        # @return [Hash] Internal features representation
        def load_features(subjectid=nil)
          if File.exist?(@uri)
            file = File.new(@uri)
          else
            file = Tempfile.new("ot-rdfxml")
            # do not concat /features to uri string, this would not work for dataset/R401577?max=3 
            uri = URI::parse(@uri)
            # PENDING
            # ambit models return http://host/dataset/id?feature_uris[]=sth but 
            # amibt dataset services does not support http://host/dataset/id/features?feature_uris[]=sth
            # and features are not inlcuded in http://host/dataset/id/features
            # -> load features from complete dataset
            uri.path = File.join(uri.path,"features") unless @uri=~/\?(feature_uris|page|pagesize)/
            uri = uri.to_s
            file.puts OpenTox::RestClientWrapper.get uri,{:subjectid => subjectid,:accept => "application/rdf+xml"},nil,false
            file.close
            to_delete = file.path
          end
          statements = []
          features = Set.new
          `rapper -i rdfxml -o ntriples #{file.path} 2>/dev/null`.each_line do |line|
            triple = line.chomp.split('> ').collect{|i| i.sub(/\s+.$/,'').gsub(/[<>"]/,'')}[0..2]
            statements << triple
            features << triple[0] if triple[1] == RDF.type and (triple[2] =~ /Feature|Substructure/) 
          end
          File.delete(to_delete) if to_delete
          statements.each do |triple|
            if features.include? triple[0]
              @dataset.features[triple[0]] = {} unless @dataset.features[triple[0]]
              if triple[1] == RDF.type
                 @dataset.features[triple[0]][triple[1]] = [] unless @dataset.features[triple[0]][triple[1]]
                 @dataset.features[triple[0]][triple[1]] << triple[2].split('^^').first
              else
                @dataset.features[triple[0]][triple[1]] = triple[2].split('^^').first
              end
            end
          end
          @dataset.features
        end

      end

    end

    # Parser for getting spreadsheet data into a dataset
    class Spreadsheets

      attr_accessor :dataset

      def initialize
        @data = []
        @features = []
        @feature_types = {}

        @format_errors = []
        @id_errors = []
        @activity_errors = []
        @duplicates = {}
        @max_class_values = 3
      end

      def detect_new_values(row, value_maps)
        row.shift
        row.each_index do |i|
          value = row[i]
          value_maps[i] = Hash.new if value_maps[i].nil?
          value_maps[i][value].nil? ? value_maps[i][value]=0 : value_maps[i][value] += 1
        end
        value_maps
      end

      # Load Spreadsheet book (created with roo gem http://roo.rubyforge.org/, excel format specification: http://toxcreate.org/help)
      # @param [Excel] book Excel workbook object (created with roo gem)
      # @return [OpenTox::Dataset] Dataset object with Excel data
      def load_spreadsheet(book, drop_missing=false)
        book.default_sheet = 0
        headers = book.row(1)
        add_features headers
        value_maps = Array.new
        regression_features=Array.new

        2.upto(book.last_row) { |i| 
          row = book.row(i)
          value_maps = detect_new_values(row, value_maps)
          value_maps.each_with_index { |vm,j|
            if vm.size > @max_class_values # 5 is the maximum nr of classes supported by Fminer.
              regression_features[j]=true 
            else
              regression_features[j]=false
            end
          }
        }

        2.upto(book.last_row) { |i| 
          drop=false
          row = book.row(i)
          raise "Entry has size #{row.size}, different from headers (#{headers.size})" if row.size != headers.size
          if row.include?("")
            @format_errors << "Row #{i} has #{row.count("")} missing values" 
            drop=true
            drop_missing=true if (row.count("") == row.size-1) 
          end
          add_values(row, regression_features) unless (drop_missing && drop)
          if (drop_missing && drop) 
            @format_errors << "Row #{i} not added" 
          end
        }
        warnings
        @dataset
      end

      # Load CSV string (format specification: http://toxcreate.org/help)
      # @param [String] csv CSV representation of the dataset
      # @return [OpenTox::Dataset] Dataset object with CSV data
      def load_csv(csv, drop_missing=false)
        row = 0
        input = csv.split("\n")
        headers = split_row(input.shift)
        add_features(headers)
        value_maps = Array.new
        regression_features=Array.new

        input.each { |row| 
          row = split_row(row)
          value_maps = detect_new_values(row, value_maps)
          value_maps.each_with_index { |vm,j|
            if vm.size > @max_class_values # max @max_class_values classes.
              regression_features[j]=true 
            else
              regression_features[j]=false
            end
          }
        }

        input.each_with_index { |row, i| 
          drop=false
          row = split_row(row)
          raise "Entry has size #{row.size}, different from headers (#{headers.size})" if row.size != headers.size
          if row.include?("")
            @format_errors << "Row #{i} has #{row.count("")} missing values" 
            drop=true
            drop_missing=true if (row.count("") == row.size-1) 
          end
          add_values(row, regression_features) unless (drop_missing && drop)
          if (drop_missing && drop) 
            @format_errors << "Row #{i} not added" 
          end
        }
        warnings
        @dataset
      end

      private

      def warnings

        info = ''
        @feature_types.each do |feature,types|
          if types.uniq.size == 0
            type = "helper#MissingFeature"
          elsif types.uniq.size > 1
            type = OT.NumericFeature
          else
            type = types.first
          end
          @dataset.add_feature_metadata(feature,{RDF.type => [type]})
          info += "\"#{@dataset.feature_name(feature)}\" detected as #{type.split('#').last}." if type

          # TODO: rewrite feature values
          # TODO if value.to_f == 0 @activity_errors << "#{id} Zero values not allowed for regression datasets - entry ignored."
        end

        @dataset.metadata[OT.Info] = info 

        warnings = ''
        warnings += "<p>Incorrect structures (ignored):</p>" + @id_errors.join("<br/>") unless @id_errors.empty?
        warnings += "<p>Irregular activities (ignored):</p>" + @activity_errors.join("<br/>") unless @activity_errors.empty?
        warnings += "<p>Format errors:</p>" + @format_errors.join("<br/>") unless @format_errors.empty?
        duplicate_warnings = ''
        @duplicates.each {|inchi,lines| duplicate_warnings << "<p>#{lines.join('<br/>')}</p>" if lines.size > 1 }
        warnings += "<p>Duplicate structures (all structures/activities used for model building, please make sure that the results were obtained from <em>independent</em> experiments):</p>" + duplicate_warnings unless duplicate_warnings.empty?

        @dataset.metadata[OT.Warnings] = warnings 

      end

      # Adds a row of features to a dataset
      # @param Array A row split up as an array
      # @return Array Indices for duplicate features
      def add_features(row)
        row=row.collect
        row.shift  # get rid of id entry
        @duplicate_feature_indices = [] # starts with 0 at first f after id
        row.each_with_index do |feature_name, idx|
          feature_uri = File.join(@dataset.uri,"feature",URI.encode(feature_name))
          unless @features.include? feature_uri
            @feature_types[feature_uri] = []
            @features << feature_uri
            @dataset.add_feature(feature_uri,{DC.title => feature_name})
          else
            @duplicate_feature_indices << idx
            @format_errors << "Duplicate Feature '#{feature_name}' at pos #{idx}"
          end
        end
      end

      # Adds a row to a dataset
      # @param Array A row split up as an array
      # @param Array Indicator for regression for each field
      # @param Array Indices for duplicate features
      def add_values(row, regression_features)

        id = row.shift
        case id
        when /InChI/
          compound = Compound.from_inchi(URI.decode_www_form_component(id))
        else
          compound = Compound.from_smiles(id)
        end

        if compound.nil? or compound.inchi.nil? or compound.inchi == ""
          @id_errors << id+", "+row.join(", ") 
          return false
        end
        @duplicates[compound.inchi] = [] unless @duplicates[compound.inchi]
        @duplicates[compound.inchi] << id+", "+row.join(", ")

        feature_idx = 0
        row.each_index do |i|

          unless @duplicate_feature_indices.include? i

            value = row[i]
            #LOGGER.warn "Missing values for #{id}" if value.size == 0 # String is empty
            feature = @features[feature_idx]
  
            type = feature_type(value) # May be NIL
            type = OT.NominalFeature unless (type.nil? || regression_features[i])
            @feature_types[feature] << type if type
  
            val = nil
            case type
            when OT.NumericFeature
              val = value.to_f
            when OT.NominalFeature
              val = value.to_s
            end

            feature_idx += 1
  
            if val != nil
              @dataset.add(compound.uri, feature, val)
              if type != OT.NumericFeature
                @dataset.features[feature][OT.acceptValue] = [] unless @dataset.features[feature][OT.acceptValue]
                @dataset.features[feature][OT.acceptValue] << val.to_s unless @dataset.features[feature][OT.acceptValue].include?(val.to_s)
              end
            end

          end

        end
      end

      def feature_type(value)
        if value == ""
          return nil
        elsif OpenTox::Algorithm::numeric? value
          return OT.NumericFeature
        else
          return OT.NominalFeature
        end
      end

      def split_row(row)
        row.chomp.gsub(/["']/,'').split(/\s*[,;\t]\s*/,-1) # -1: do not skip empty cells
      end

    end

    class Table

      attr_accessor :data, :features, :compounds

      def initialize
        @data = {}
        @activity_errors = []
        @max_class_values = 3
      end

      def feature_values(feature)
        @data.collect{|c, row| row[feature]}.uniq.compact
      end

      def feature_types(feature)
        @data.collect{|c, row| feature_type(row[feature])}.uniq.compact
      end

      def features
        @data.collect{|c,row| row.keys}.flatten.uniq
      end

      def clean_features
        ignored_features = []
        features.each do |feature|
          if feature_values(feature).size > @max_class_values
            if feature_types(feature).size == 1 and feature_types(feature).first == OT.NumericFeature
              # REGRESSION
            elsif feature_types(feature).include? OT.NumericFeature
              @data.each{|c,row| row[feature] = nil unless OpenTox::Algorithm::numeric?(row[feature]) } # delete nominal features
              @activity_errors << "Nominal feature values of #{feature} ignored (using numeric features for regression models)."
            else
              @activity_errors << "Feature #{feature} ignored (more than #{@max_class_values} nominal feature values and no numeric values)."
              ignored_features << feature
              next
            end
          elsif feature_values(feature).size <= 1
              @activity_errors << "Feature #{feature} ignored (less than 2 feature values)."
              ignored_features << feature
          else
            # CLASSIFICATION
          end
        end
        ignored_features.each do |feature|
          @data.each{ |c,row| row.delete feature }
        end
        @activity_errors
      end

      def add_to_dataset(dataset)
        features.each do |feature_name|
          feature_uri = File.join(dataset.uri,"feature",URI.encode(feature_name))
          dataset.add_feature(feature_uri,{DC.title => feature_name})
        end

        @data.each do |compound,row|
          unless row.empty?
            row.each do |feature,value|
              if OpenTox::Algorithm::numeric?(value)
                value = value.to_f
              elsif value.nil? or value.empty?
                value = nil
              else
                value = value.to_s
              end
              feature_uri = File.join(dataset.uri,"feature",URI.encode(feature))
              dataset.add(compound, feature_uri, value)
              #dataset.features[feature_uri][RDF.type] = feature_types(feature)
              #dataset.features[feature_uri][OT.acceptValue] = feature_values(feature)
              if feature_types(feature).include? OT.NumericFeature
                dataset.features[feature_uri][RDF.type] = [OT.NumericFeature]
              else
                dataset.features[feature_uri][RDF.type] = [OT.NominalFeature]
                dataset.features[feature_uri][OT.acceptValue] = feature_values(feature) 
              end
            end
          end
        end
      end

      private

      def feature_type(value)
        if value.nil?
          return nil
        elsif OpenTox::Algorithm::numeric? value
          return OT.NumericFeature
        else
          return OT.NominalFeature
        end
      end

    end

    # quick hack to enable sdf import via csv
    # should be refactored 
    class Sdf

      attr_accessor :dataset

      def initialize
        @data = {}

        @compound_errors = []
        @activity_errors = []
        @duplicates = {}
      end

      def load_sdf(sdf)

        obconversion = OpenBabel::OBConversion.new
        obmol = OpenBabel::OBMol.new
        obconversion.set_in_and_out_formats "sdf", "inchi"

        table = Table.new

        properties = []
        sdf.each_line { |l| properties << l.to_s if l.match(/</) }
        properties.uniq!
        properties.sort!
        properties.collect!{ |p| p.gsub(/<|>/,'').strip.chomp }

        rec = 0
        sdf.split(/\$\$\$\$\r*\n/).each do |s|
          rec += 1
          obconversion.read_string obmol, s
          begin
            inchi = obconversion.write_string(obmol).gsub(/\s/,'').chomp 
            @duplicates[inchi] = [] unless @duplicates[inchi]
            @duplicates[inchi] << rec #inchi#+", "+row.join(", ")
            compound = Compound.from_inchi inchi
          rescue
            @compound_errors << "Could not convert structure to InChI, all entries for this compound (record #{rec}) have been ignored! \n#{s}"
            next
          end
          row = {}
          obmol.get_data.each { |d| row[d.get_attribute] = d.get_value if properties.include?(d.get_attribute) }
          table.data[compound.uri] = row
        end
        
        # find and remove ignored_features
        @activity_errors = table.clean_features
        table.add_to_dataset @dataset

        warnings = ''
        warnings += "<p>Incorrect structures (ignored):</p>" + @compound_errors.join("<br/>") unless @compound_errors.empty?
        warnings += "<p>Irregular activities (ignored):</p>" + @activity_errors.join("<br/>") unless @activity_errors.empty?
        duplicate_warnings = ''
        @duplicates.each {|inchi,lines| duplicate_warnings << "<p>#{lines.join('<br/>')}</p>" if lines.size > 1 }
        warnings += "<p>Duplicated structures (all structures/activities used for model building, please  make sure, that the results were obtained from <em>independent</em> experiments):</p>" + duplicate_warnings unless duplicate_warnings.empty?

        @dataset.metadata[OT.Warnings] = warnings 
        @dataset

      end

    end
  end
end
