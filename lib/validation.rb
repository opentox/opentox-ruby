require "yaml"
module OpenTox
  class Validation
    include OpenTox
    
    # find validation, raises error if not found
    # @param [String] uri
    # @param [String,optional] subjectid
    # @return [OpenTox::Validation]
    def self.find( uri, subjectid=nil )
      val = Validation.new(uri)
      val.load_metadata( subjectid )
      val
    end
    
    # returns a filtered list of validation uris
    # @param [Hash,optional] params, validation-params to filter the uris (could be model, training_dataset, ..)  
    # @return [Array]    
    def self.list( params={} )
      filter_string = ""
      params.each do |k,v|
        filter_string = "?" if filter_string.length==0 
        filter_string += k.to_s+"="+v
      end
      (OpenTox::RestClientWrapper.get(CONFIG[:services]["opentox-validation"]+filter_string).split("\n"))
    end
    
    # creates a training test split validation, waits until it finishes, may take some time
    # @param [Hash] params (required:algorithm_uri,dataset_uri,prediction_feature, optional:algorithm_params,split_ratio(0.67),random_seed(1))
    # @param [String,optional] subjectid
    # @param [OpenTox::Task,optional] waiting_task (can be a OpenTox::Subtask as well), progress is updated accordingly
    # @return [OpenTox::Validation]
    def self.create_training_test_split( params, subjectid=nil, waiting_task=nil )
      params[:subjectid] = subjectid if subjectid
      uri = OpenTox::RestClientWrapper.post( File.join(CONFIG[:services]["opentox-validation"],"training_test_split"),
        params,{:content_type => "text/uri-list"},waiting_task )
      Validation.new(uri)
    end
    
    # creates a training test validation, waits until it finishes, may take some time
    # @param [Hash] params (required:algorithm_uri,training_dataset_uri,prediction_feature,test_dataset_uri,optional:algorithm_params)
    # @param [String,optional] subjectid
    # @param [OpenTox::Task,optional] waiting_task (can be a OpenTox::Subtask as well), progress is updated accordingly
    # @return [OpenTox::Validation]
    def self.create_training_test_validation( params, subjectid=nil, waiting_task=nil )
      params[:subjectid] = subjectid if subjectid
      uri = OpenTox::RestClientWrapper.post( File.join(CONFIG[:services]["opentox-validation"],"training_test_validation"),
        params,{:content_type => "text/uri-list"},waiting_task )
      Validation.new(uri)
    end
    
    # creates a bootstrapping validation, waits until it finishes, may take some time
    # @param [Hash] params (required:algorithm_uri,dataset_uri,prediction_feature, optional:algorithm_params,random_seed(1))
    # @param [String,optional] subjectid
    # @param [OpenTox::Task,optional] waiting_task (can be a OpenTox::Subtask as well), progress is updated accordingly
    # @return [OpenTox::Validation]
    def self.create_bootstrapping_validation( params, subjectid=nil, waiting_task=nil )
      params[:subjectid] = subjectid if subjectid
      uri = OpenTox::RestClientWrapper.post( File.join(CONFIG[:services]["opentox-validation"],"bootstrapping"),
        params,{:content_type => "text/uri-list"},waiting_task )
      Validation.new(uri)
    end
    
    # looks for report for this validation, creates a report if no report is found
    # @param [String,optional] subjectid
    # @param [OpenTox::Task,optional] waiting_task (can be a OpenTox::Subtask as well), progress is updated accordingly
    # @return [String] report uri
    def find_or_create_report( subjectid=nil, waiting_task=nil )
      @report = ValidationReport.find_for_validation(@uri, subjectid) unless @report
      @report = ValidationReport.create(@uri, {}, subjectid, waiting_task) unless @report
      @report.uri
    end
    
    # creates a validation object from crossvaldiation statistics, raise error if not found
    # (as crossvaldiation statistics are returned as an average valdidation over all folds)
    # @param [String] crossvalidation uri
    # @param [String,optional] subjectid
    # @return [OpenTox::Validation]
    def self.from_cv_statistics( crossvalidation_uri, subjectid=nil )
      find( File.join(crossvalidation_uri, 'statistics'),subjectid )
    end
    
    # loads metadata via yaml from validation object
    # fields (like for example the validated model) can be acces via validation.metadata[OT.model]
    def load_metadata( subjectid=nil )
      @metadata = YAML.load(OpenTox::RestClientWrapper.get(uri,{:subjectid => subjectid, :accept => "application/x-yaml"}))
    end
    
    # returns confusion matrix as array, predicted values are in rows
    # example:
    # [[nil,"active","moderate","inactive"],["active",1,3,99],["moderate",4,2,8],["inactive",3,8,6]]
    # -> 99 inactive compounds have been predicted as active 
    def confusion_matrix
      raise "no classification statistics, probably a regression valdiation" unless @metadata[OT.classificationStatistics]
      matrix =  @metadata[OT.classificationStatistics][OT.confusionMatrix][OT.confusionMatrixCell]
      values = matrix.collect{|cell| cell[OT.confusionMatrixPredicted]}.uniq
      table = [[nil]+values]
      values.each do |c|
        table << [c]
        values.each do |r|
          matrix.each do |cell|
            if cell[OT.confusionMatrixPredicted]==c and cell[OT.confusionMatrixActual]==r
              table[-1] << cell[OT.confusionMatrixValue].to_f
              break
            end
          end
        end
      end
      table
    end
    
    # returns probability-distribution for a given prediction
    # it takes all predictions into account that have a confidence value that is >= confidence and that have the same predicted value
    # (minimum 12 predictions with the hightest confidence are selected (even if the confidence is lower than the given param)
    # 
    # @param [Float] confidence value (between 0 and 1)
    # @param [String] predicted value
    # @param [String,optional] subjectid
    # @return [Hash] see example
    #
    # Example 1:
    # validation.probabilities(0.3,"active")
    # -> {:min_confidence=>0.32, :num_predictions=>20, :probs=>{"active"=>0.7, "moderate"=>0.25 "inactive"=>0.05}}
    # there have been 20 "active" predictions with confidence >= 0.3, 70 percent of them beeing correct
    #
    # Example 2:
    # validation.probabilities(0.8,"active")
    # -> {:min_confidence=>0.45, :num_predictions=>12, :probs=>{"active"=>0.9, "moderate"=>0.1 "inactive"=>0}}
    # the given confidence value was to high (i.e. <12 predictions with confidence value >= 0.8)
    # the top 12 "active" predictions have a min_confidence of 0.45, 90 percent of them beeing correct
    # 
    def probabilities( confidence, prediction, subjectid=nil )
      YAML.load(OpenTox::RestClientWrapper.get(@uri+"/probabilities?prediction="+prediction.to_s+"&confidence="+confidence.to_s,
        {:subjectid => subjectid, :accept => "application/x-yaml"}))
    end
  end
  
  class Crossvalidation
    include OpenTox

    attr_reader :report
    
    # find crossvalidation, raises error if not found
    # @param [String] uri
    # @param [String,optional] subjectid
    # @return [OpenTox::Crossvalidation]
    def self.find( uri, subjectid=nil )
      cv = Crossvalidation.new(uri)
      cv.load_metadata( subjectid )
      cv
    end
    
    # returns a filtered list of crossvalidation uris
    # @param [Hash,optional] params, crossvalidation-params to filter the uris (could be algorithm, dataset, ..)  
    # @return [Array]    
    def self.list( params={} )
      filter_string = ""
      params.each do |k,v|
        filter_string = "?" if filter_string.length==0 
        filter_string += k.to_s+"="+v
      end
      (OpenTox::RestClientWrapper.get(File.join(CONFIG[:services]["opentox-validation"],"crossvalidation")+filter_string).split("\n"))
    end
		
    # creates a crossvalidations, waits until it finishes, may take some time
    # @param [Hash] params (required:algorithm_uri,dataset_uri,prediction_feature, optional:algorithm_params,num_folds(10),random_seed(1),stratified(false))
    # @param [String,optional] subjectid
    # @param [OpenTox::Task,optional] waiting_task (can be a OpenTox::Subtask as well), progress is updated accordingly
    # @return [OpenTox::Crossvalidation]
    def self.create( params, subjectid=nil, waiting_task=nil )
      params[:subjectid] = subjectid if subjectid
      uri = OpenTox::RestClientWrapper.post( File.join(CONFIG[:services]["opentox-validation"],"crossvalidation"),
        params,{:content_type => "text/uri-list"},waiting_task )
      Crossvalidation.new(uri)
    end

    # looks for report for this crossvalidation, creates a report if no report is found
    # @param [String,optional] subjectid
    # @param [OpenTox::Task,optional] waiting_task (can be a OpenTox::Subtask as well), progress is updated accordingly
    # @return [String] report uri
    def find_or_create_report( subjectid=nil, waiting_task=nil )
      @report = CrossvalidationReport.find_for_crossvalidation(@uri, subjectid) unless @report
      @report = CrossvalidationReport.create(@uri, subjectid, waiting_task) unless @report
      @report.uri
    end
    
    # loads metadata via yaml from crossvalidation object
    # fields (like for example the validations) can be acces via validation.metadata[OT.validation]
    def load_metadata( subjectid=nil )
      @metadata = YAML.load(OpenTox::RestClientWrapper.get(uri,{:subjectid => subjectid, :accept => "application/x-yaml"}))
    end
    
    # returns a Validation object containing the statistics of the crossavlidation
    def statistics( subjectid=nil )
      Validation.from_cv_statistics( @uri, subjectid )
    end
    
    # documentation see OpenTox::Validation.probabilities
    def probabilities( confidence, prediction, subjectid=nil )
      YAML.load(OpenTox::RestClientWrapper.get(@uri+"/statistics/probabilities?prediction="+prediction.to_s+"&confidence="+confidence.to_s,
        {:subjectid => subjectid, :accept => "application/x-yaml"}))
    end
    
  end
  
  class ValidationReport
    include OpenTox
    
    # finds ValidationReport via uri, raises error if not found
    # @param [String] uri
    # @param [String,optional] subjectid
    # @return [OpenTox::ValidationReport]
    def self.find( uri, subjectid=nil )
      OpenTox::RestClientWrapper.get(uri,{:subjectid => subjectid})
      rep = ValidationReport.new(uri)
      rep.load_metadata( subjectid )
      rep
    end
    
    # finds ValidationReport for a particular validation
    # @param [String] crossvalidation uri 
    # @param [String,optional] subjectid
    # @return [OpenTox::ValidationReport] nil if no report found
    def self.find_for_validation( validation_uri, subjectid=nil )
      uris = RestClientWrapper.get(File.join(CONFIG[:services]["opentox-validation"],
        "/report/validation?validation="+validation_uri), {:subjectid => subjectid}).chomp.split("\n")
      uris.size==0 ? nil : ValidationReport.new(uris[-1])
    end
    
    # creates a validation report via validation
    # @param [String] validation uri 
    # @param [Hash] params addiditonal possible 
    #               (min_confidence, params={}, min_num_predictions, max_num_predictions)
    # @param [String,optional] subjectid
    # @param [OpenTox::Task,optional] waiting_task (can be a OpenTox::Subtask as well), progress is updated accordingly
    # @return [OpenTox::ValidationReport]
    def self.create( validation_uri, params={}, subjectid=nil, waiting_task=nil )
      params = {} if params==nil
      raise OpenTox::BadRequestError.new "params is no hash" unless params.is_a?(Hash)
      params[:validation_uris] = validation_uri
      params[:subjectid] = subjectid
      uri = RestClientWrapper.post(File.join(CONFIG[:services]["opentox-validation"],"/report/validation"),
        params, {}, waiting_task )
      ValidationReport.new(uri)
    end
    
  end

  class CrossvalidationReport
    include OpenTox
    
    # finds CrossvalidationReport via uri, raises error if not found
    # @param [String] uri
    # @param [String,optional] subjectid
    # @return [OpenTox::CrossvalidationReport]
    def self.find( uri, subjectid=nil )
      OpenTox::RestClientWrapper.get(uri,{:subjectid => subjectid})
      rep = CrossvalidationReport.new(uri)
      rep.load_metadata( subjectid )
      rep
    end
    
    # finds CrossvalidationReport for a particular crossvalidation
    # @param [String] crossvalidation uri 
    # @param [String,optional] subjectid
    # @return [OpenTox::CrossvalidationReport] nil if no report found
    def self.find_for_crossvalidation( crossvalidation_uri, subjectid=nil )
      uris = RestClientWrapper.get(File.join(CONFIG[:services]["opentox-validation"],
        "/report/crossvalidation?crossvalidation="+crossvalidation_uri), {:subjectid => subjectid}).chomp.split("\n")
      uris.size==0 ? nil : CrossvalidationReport.new(uris[-1])
    end
    
    # creates a crossvalidation report via crossvalidation
    # @param [String] crossvalidation uri 
    # @param [String,optional] subjectid
    # @param [OpenTox::Task,optional] waiting_task (can be a OpenTox::Subtask as well), progress is updated accordingly
    # @return [OpenTox::CrossvalidationReport]
    def self.create( crossvalidation_uri, subjectid=nil, waiting_task=nil )
      uri = RestClientWrapper.post(File.join(CONFIG[:services]["opentox-validation"],"/report/crossvalidation"),
        { :validation_uris => crossvalidation_uri, :subjectid => subjectid }, {}, waiting_task )
      CrossvalidationReport.new(uri)
    end
  end
  
  
  class AlgorithmComparisonReport
    include OpenTox
    
    # finds AlgorithmComparisonReport via uri, raises error if not found
    # @param [String] uri
    # @param [String,optional] subjectid
    # @return [OpenTox::CrossvalidationReport]
    def self.find( uri, subjectid=nil )
      OpenTox::RestClientWrapper.get(uri,{:subjectid => subjectid})
      rep = AlgorithmComparisonReport.new(uri)
      rep.load_metadata( subjectid )
      rep
    end
    
    # finds AlgorithmComparisonReport for a particular crossvalidation
    # @param [String] crossvalidation uri 
    # @param [String,optional] subjectid
    # @return [OpenTox::AlgorithmComparisonReport] nil if no report found
    def self.find_for_crossvalidation( crossvalidation_uri, subjectid=nil )
      uris = RestClientWrapper.get(File.join(CONFIG[:services]["opentox-validation"],
        "/report/algorithm_comparison?crossvalidation="+crossvalidation_uri), {:subjectid => subjectid}).chomp.split("\n")
      uris.size==0 ? nil : AlgorithmComparisonReport.new(uris[-1])
    end
    
    # creates a algorithm comparison report via crossvalidation uris
    # @param [Hash] crossvalidation uri_hash, see example 
    # @param [Hash] params addiditonal possible 
    #               (ttest_significance, ttest_attributes, min_confidence, min_num_predictions, max_num_predictions)
    # @param [String,optional] subjectid
    # @param [OpenTox::Task,optional] waiting_task (can be a OpenTox::Subtask as well), progress is updated accordingly
    # @return [OpenTox::AlgorithmComparisonReport]
    # example for hash:
    # { :lazar-bbrc => [ http://host/validation/crossvalidation/x1, http://host/validation/crossvalidation/x2 ],
    #   :lazar-last => [ http://host/validation/crossvalidation/xy, http://host/validation/crossvalidation/xy ] }
    def self.create( crossvalidation_uri_hash, params={}, subjectid=nil, waiting_task=nil )
      identifier = []
      validation_uris = []
      crossvalidation_uri_hash.each do |id, uris|
        uris.each do |uri|
          identifier << id
          validation_uris << uri
        end
      end
      params = {} if params==nil
      raise OpenTox::BadRequestError.new "params is no hash" unless params.is_a?(Hash)
      params[:validation_uris] = validation_uris.join(",")
      params[:identifier] = identifier.join(",")
      params[:subjectid] = subjectid
      uri = RestClientWrapper.post(File.join(CONFIG[:services]["opentox-validation"],"/report/algorithm_comparison"),
        params, {}, waiting_task )
      AlgorithmComparisonReport.new(uri)
    end
  end  
  
  class QMRFReport
    include OpenTox
    
    # finds QMRFReport, raises Error if not found
    # @param [String] uri
    # @param [String,optional] subjectid
    # @return [OpenTox::QMRFReport]
    def self.find( uri, subjectid=nil )
      # PENDING load crossvalidation data?
      OpenTox::RestClientWrapper.get(uri,{:subjectid => subjectid})
      QMRFReport.new(uri)
    end
    
    # finds QMRF report for a particular model
    # @param [String] model_uri 
    # @param [String,optional] subjectid
    # @return [OpenTox::QMRFReport] nil if no report found
    def self.find_for_model( model_uri, subjectid=nil )
      uris = RestClientWrapper.get(File.join(CONFIG[:services]["opentox-validation"],
        "/reach_report/qmrf?model="+model_uri), {:subjectid => subjectid}).chomp.split("\n")
      uris.size==0 ? nil : QMRFReport.new(uris[-1])
    end
    
    # creates a qmrf report via model
    # @param [String] model_uri 
    # @param [String,optional] subjectid
    # @param [OpenTox::Task,optional] waiting_task (can be a OpenTox::Subtask as well), progress is updated accordingly
    # @return [OpenTox::QMRFReport]
    def self.create( model_uri, subjectid=nil, waiting_task=nil )
      uri = RestClientWrapper.post(File.join(CONFIG[:services]["opentox-validation"],"/reach_report/qmrf"), 
        { :model_uri => model_uri, :subjectid => subjectid }, {}, waiting_task )
      QMRFReport.new(uri)
    end
  end
  
end

