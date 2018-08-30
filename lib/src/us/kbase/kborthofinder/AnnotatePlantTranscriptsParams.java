
package us.kbase.kborthofinder;

import java.util.HashMap;
import java.util.Map;
import javax.annotation.Generated;
import com.fasterxml.jackson.annotation.JsonAnyGetter;
import com.fasterxml.jackson.annotation.JsonAnySetter;
import com.fasterxml.jackson.annotation.JsonInclude;
import com.fasterxml.jackson.annotation.JsonProperty;
import com.fasterxml.jackson.annotation.JsonPropertyOrder;


/**
 * <p>Original spec-file type: AnnotatePlantTranscriptsParams</p>
 * 
 * 
 */
@JsonInclude(JsonInclude.Include.NON_NULL)
@Generated("com.googlecode.jsonschema2pojo")
@JsonPropertyOrder({
    "threshold",
    "input_ws",
    "input_genome",
    "output_genome"
})
public class AnnotatePlantTranscriptsParams {

    @JsonProperty("threshold")
    private Double threshold;
    @JsonProperty("input_ws")
    private String inputWs;
    @JsonProperty("input_genome")
    private String inputGenome;
    @JsonProperty("output_genome")
    private String outputGenome;
    private Map<String, Object> additionalProperties = new HashMap<String, Object>();

    @JsonProperty("threshold")
    public Double getThreshold() {
        return threshold;
    }

    @JsonProperty("threshold")
    public void setThreshold(Double threshold) {
        this.threshold = threshold;
    }

    public AnnotatePlantTranscriptsParams withThreshold(Double threshold) {
        this.threshold = threshold;
        return this;
    }

    @JsonProperty("input_ws")
    public String getInputWs() {
        return inputWs;
    }

    @JsonProperty("input_ws")
    public void setInputWs(String inputWs) {
        this.inputWs = inputWs;
    }

    public AnnotatePlantTranscriptsParams withInputWs(String inputWs) {
        this.inputWs = inputWs;
        return this;
    }

    @JsonProperty("input_genome")
    public String getInputGenome() {
        return inputGenome;
    }

    @JsonProperty("input_genome")
    public void setInputGenome(String inputGenome) {
        this.inputGenome = inputGenome;
    }

    public AnnotatePlantTranscriptsParams withInputGenome(String inputGenome) {
        this.inputGenome = inputGenome;
        return this;
    }

    @JsonProperty("output_genome")
    public String getOutputGenome() {
        return outputGenome;
    }

    @JsonProperty("output_genome")
    public void setOutputGenome(String outputGenome) {
        this.outputGenome = outputGenome;
    }

    public AnnotatePlantTranscriptsParams withOutputGenome(String outputGenome) {
        this.outputGenome = outputGenome;
        return this;
    }

    @JsonAnyGetter
    public Map<String, Object> getAdditionalProperties() {
        return this.additionalProperties;
    }

    @JsonAnySetter
    public void setAdditionalProperties(String name, Object value) {
        this.additionalProperties.put(name, value);
    }

    @Override
    public String toString() {
        return ((((((((((("AnnotatePlantTranscriptsParams"+" [threshold=")+ threshold)+", inputWs=")+ inputWs)+", inputGenome=")+ inputGenome)+", outputGenome=")+ outputGenome)+", additionalProperties=")+ additionalProperties)+"]");
    }

}
