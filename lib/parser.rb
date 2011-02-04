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
        if @dataset
          uri = File.join(@uri,"metadata")
        else
          uri = @uri
        end
        uri += "?subjectid=#{CGI.escape(subjectid)}" if subjectid 
        statements = []
        parameter_ids = []
        `rapper -i rdfxml -o ntriples #{uri} 2>/dev/null`.each_line do |line|
          triple = line.to_triple
          @metadata[triple[1]] = triple[2].split('^^').first if triple[0] == @uri and triple[1] != RDF['type']
          statements << triple 
          parameter_ids << triple[2] if triple[1] == OT.parameters
        end
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
      
      # loads metadata from rdf-data
      # @param [String] rdf
      # @param [String] type of the info (e.g. OT.Task, OT.ErrorReport) needed to get the subject-uri
      # @return [Hash] metadata 
      def self.metadata_from_rdf( rdf, type )
        # write to file and read convert with rapper into tripples
        file = Tempfile.new("ot-rdfxml")
        file.puts rdf
        file.close
        file = "file://"+file.path
        #puts "cmd: rapper -i rdfxml -o ntriples #{file} 2>/dev/null"
        triples = `rapper -i rdfxml -o ntriples #{file} 2>/dev/null`
        
        # load uri via type
        uri = nil
        triples.each_line do |line|
          triple = line.to_triple
          if triple[1] == RDF['type'] and triple[2]==type
             raise "uri already set, two uris found with type: "+type.to_s if uri
             uri = triple[0]
          end
        end
        
        # load metadata
        metadata = {}
        triples.each_line do |line|
          triple = line.to_triple
          metadata[triple[1]] = triple[2].split('^^').first if triple[0] == uri and triple[1] != RDF['type']
        end
        metadata
      end

      # Generic parser for all OpenTox classes
      class Generic
        include Owl
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
          uri = @uri
          uri += "?subjectid=#{CGI.escape(subjectid)}" if subjectid
          data = {}
          feature_values = {}
          feature = {}
          other_statements = {}
          `rapper -i rdfxml -o ntriples #{uri} 2>/dev/null`.each_line do |line|
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
            else 
            end
          end
          data.each do |id,entry|
            entry[:values].each do |value_id|
              value = feature_values[value_id].split(/\^\^/).first # remove XSD.type
              @dataset.add entry[:compound],feature[value_id],value
            end
          end
          load_features
          @dataset.metadata = load_metadata
          @dataset
        end

        # Read only features from a dataset service. 
        # @return [Hash] Internal features representation
        def load_features(subjectid=nil)
          uri = File.join(@uri,"features")
          uri += "?subjectid=#{CGI.escape(subjectid)}" if subjectid 
          statements = []
          features = Set.new
          `rapper -i rdfxml -o ntriples #{uri} 2>/dev/null`.each_line do |line|
            triple = line.chomp.split('> ').collect{|i| i.sub(/\s+.$/,'').gsub(/[<>"]/,'')}[0..2]
            statements << triple
            features << triple[0] if triple[1] == RDF['type'] and triple[2] == OT.Feature
          end
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
        2.upto(book.last_row) { |i| add_values book.row(i) }
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
        input.each { |row| add_values split_row(row) }
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
          @dataset.add_feature_metadata(feature,{OT.isA => type})
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

      def add_values(row)

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

          case type
          when OT.NominalFeature
            case value.to_s
            when TRUE_REGEXP
              @dataset.add(compound.uri, feature, true )
            when FALSE_REGEXP
              @dataset.add(compound.uri, feature, false )
            end
          when OT.NumericFeature
            @dataset.add compound.uri, feature, value.to_f
          when OT.StringFeature
            @dataset.add compound.uri, feature, value.to_s
            @activity_errors << smiles+", "+row.join(", ")
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
