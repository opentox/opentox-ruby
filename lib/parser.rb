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
            # do not concat /metadata to uri string, this would not work for dataset/R401577?max=3 
            uri = URI::parse(@uri)
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
            if triple[1] == RDF.type # allow multiple types
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
        @metadata
      end
      
      # creates owl object from rdf-data
      # @param [String] rdf
      # @param [String] type of the info (e.g. OT.Task, OT.ErrorReport) needed to get the subject-uri
      # @return [Owl] with uri and metadata set 
      def self.from_rdf( rdf, type )
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
             raise "uri already set, two uris found with type: "+type.to_s if uri
             uri = triple[0]
          end
        end
        File.delete(file.path)
        # load metadata
        metadata = {}
        triples.each_line do |line|
          triple = line.to_triple
          metadata[triple[1]] = triple[2].split('^^').first if triple[0] == uri and triple[1] != RDF['type']
        end
        owl = Owl::Generic.new(uri)
        owl.metadata = metadata
        owl
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
            uri.path = File.join(uri.path,"features")
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
              @dataset.features[triple[0]][triple[1]] = triple[2].split('^^').first
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

        @format_errors = ""
        @smiles_errors = []
        @activity_errors = []
        @duplicates = {}
      end

      # Load Spreadsheet book (created with roo gem http://roo.rubyforge.org/, excel format specification: http://toxcreate.org/help)
      # @param [Excel] book Excel workbook object (created with roo gem)
      # @return [OpenTox::Dataset] Dataset object with Excel data
      def load_spreadsheet(book)
        book.default_sheet = 0
        add_features book.row(1)

        # AM: fix mixed read in
        regression_features=false
        2.upto(book.last_row) { |i| 
          row = book.row(i)
          regression_features = detect_regression_features row
          break if regression_features==true
        }
        
        2.upto(book.last_row) { |i| add_values book.row(i),regression_features }
        warnings
        @dataset
      end

      # Load CSV string (format specification: http://toxcreate.org/help)
      # @param [String] csv CSV representation of the dataset
      # @return [OpenTox::Dataset] Dataset object with CSV data
      def load_csv(csv)
        row = 0
        input = csv.split("\n")
        add_features split_row(input.shift)


        # AM: fix mixed read in
        regression_features=false
        input.each { |row| 
          row = split_row(row)
          regression_features = detect_regression_features row
          break if regression_features==true
        }
        input.each { |row| add_values split_row(row),regression_features }
        warnings
        @dataset
      end


      private

      def warnings

        info = ''
        @feature_types.each do |feature,types|
          if types.uniq.size > 1
            type = OT.NumericFeature
          else
            type = types.first
          end
          @dataset.add_feature_metadata(feature,{RDF.type => [type]})
          info += "\"#{@dataset.feature_name(feature)}\" detected as #{type.split('#').last}."

          # TODO: rewrite feature values
          # TODO if value.to_f == 0 @activity_errors << "#{smiles} Zero values not allowed for regression datasets - entry ignored."
        end

        @dataset.metadata[OT.Info] = info 

        warnings = ''
        warnings += "<p>Incorrect Smiles structures (ignored):</p>" + @smiles_errors.join("<br/>") unless @smiles_errors.empty?
        warnings += "<p>Irregular activities (ignored):</p>" + @activity_errors.join("<br/>") unless @activity_errors.empty?
        duplicate_warnings = ''
        @duplicates.each {|inchi,lines| duplicate_warnings << "<p>#{lines.join('<br/>')}</p>" if lines.size > 1 }
        warnings += "<p>Duplicated structures (all structures/activities used for model building, please  make sure, that the results were obtained from <em>independent</em> experiments):</p>" + duplicate_warnings unless duplicate_warnings.empty?

        @dataset.metadata[OT.Warnings] = warnings 

      end

      def add_features(row)
        row.shift  # get rid of smiles entry
        row.each do |feature_name|
          feature_uri = File.join(@dataset.uri,"feature",URI.encode(feature_name))
          @feature_types[feature_uri] = []
          @features << feature_uri
          @dataset.add_feature(feature_uri,{DC.title => feature_name})
        end
      end

      def detect_regression_features row
        row.shift
        regression_features=false
        row.each_index do |i|
          value = row[i]
          type = feature_type(value)
          if type == OT.NumericFeature
            regression_features=true
          end
        end
        regression_features
      end

      def add_values(row, regression_features=false)

        smiles = row.shift
        compound = Compound.from_smiles(smiles)
        if compound.nil? or compound.inchi.nil? or compound.inchi == ""
          @smiles_errors << smiles+", "+row.join(", ") 
          return false
        end
        @duplicates[compound.inchi] = [] unless @duplicates[compound.inchi]
        @duplicates[compound.inchi] << smiles+", "+row.join(", ")

        row.each_index do |i|
          value = row[i]
          feature = @features[i]
          type = feature_type(value)

          @feature_types[feature] << type 

          if (regression_features)
            val = value.to_f
          else
            case type
            when OT.NominalFeature
              case value.to_s
              when TRUE_REGEXP
                val = true
              when FALSE_REGEXP
                val = false
              end
            when OT.NumericFeature
              val = value.to_f
            when OT.StringFeature
              val = value.to_s
              @activity_errors << smiles+", "+row.join(", ")
            end
          end
          if val!=nil
            @dataset.add(compound.uri, feature, val)
            if type!=OT.NumericFeature
              @dataset.features[feature][OT.acceptValue] = [] unless @dataset.features[feature][OT.acceptValue]
              @dataset.features[feature][OT.acceptValue] << val.to_s unless @dataset.features[feature][OT.acceptValue].include?(val.to_s)
            end
          end
        end
      end

      def numeric?(value)
        true if Float(value) rescue false
      end

      def classification?(value)
        !value.to_s.strip.match(TRUE_REGEXP).nil? or !value.to_s.strip.match(FALSE_REGEXP).nil?
      end

      def feature_type(value)
        if classification? value
          return OT.NominalFeature
        elsif numeric? value
          return OT.NumericFeature
        else
          return OT.StringFeature
        end
      end

      def split_row(row)
        row.chomp.gsub(/["']/,'').split(/\s*[,;]\s*/) # remove quotes
      end

    end
  end
end
