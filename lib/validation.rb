module OpenTox
	class Validation
    include OpenTox

		attr_accessor :report_uri, :qmrf_report_uri
		
		def self.create_crossvalidation(params)
			params[:uri] = File.join(CONFIG[:services]['opentox-validation'], "crossvalidation")
      params[:num_folds] = 10 unless params[:num_folds]
		 	params[:random_seed] = 2 unless params[:random_seed]
		 	params[:stratified] = false unless params[:stratified]
      uri = OpenTox::RestClientWrapper.post(File.join(CONFIG[:services]["opentox-validation"],"/crossvalidation"),params,nil,false)
			OpenTox::Validation.new(uri)
		end

    def create_report(subjectid=nil)
			@report_uri = OpenTox::RestClientWrapper.post(File.join(CONFIG[:services]["opentox-validation"],"/report/crossvalidation"), {:validation_uris => @uri, :subjectid => subjectid}).to_s
      @report_uri
    end

    def create_qmrf_report(subjectid=nil)
      @qmrf_report_uri = OpenTox::RestClientWrapper.post(File.join(CONFIG[:services]["opentox-validation"],"/reach_report/qmrf"), {:model_uri => @uri, :subjectid => subjectid}).to_s
      @qmrf_report_uri
    end

    def summary(type, subjectid=nil)
      v = YAML.load OpenTox::RestClientWrapper.get(File.join(@uri, 'statistics'),{:accept => "application/x-yaml", :subjectid => subjectid}).to_s

      case type
      when "classification"
        tp=0; tn=0; fp=0; fn=0; n=0
        v[:classification_statistics][:confusion_matrix][:confusion_matrix_cell].each do |cell|
          if cell[:confusion_matrix_predicted] == "true" and cell[:confusion_matrix_actual] == "true"
            tp = cell[:confusion_matrix_value]
            n += tp
          elsif cell[:confusion_matrix_predicted] == "false" and cell[:confusion_matrix_actual] == "false"
            tn = cell[:confusion_matrix_value]
            n += tn
          elsif cell[:confusion_matrix_predicted] == "false" and cell[:confusion_matrix_actual] == "true"
            fn = cell[:confusion_matrix_value]
            n += fn
          elsif cell[:confusion_matrix_predicted] == "true" and cell[:confusion_matrix_actual] == "false"
            fp = cell[:confusion_matrix_value]
            n += fp
          end
        end
        {
          :nr_predictions => n,
          :true_positives => tp,
          :false_positives => fp,
          :true_negatives => tn,
          :false_negatives => fn,
          :correct_predictions => 100*(tp+tn).to_f/n,
          :weighted_area_under_roc => v[:classification_statistics][:weighted_area_under_roc].to_f,
          :sensitivity => tp.to_f/(tp+fn),
          :specificity => tn.to_f/(tn+fp),
        }
      when "regression"
        {
          :nr_predictions => v[:num_instances] - v[:num_unpredicted],
          :r_square => v[:regression_statistics][:r_square],
          :root_mean_squared_error => v[:regression_statistics][:root_mean_squared_error],
          :mean_absolute_error => v[:regression_statistics][:mean_absolute_error],
        }
      end
    end

	end
end

