import java.util.*;

public class Gene {

	public int geneId;
	public List<Double> expression;

	public int getGeneId() {
		return geneId;
	}

	public void setGeneId(int geneId) {
		this.geneId = geneId;
	}

	public List<Double> getExpression() {
		return expression;
	}

	public void setExpression(List<Double> expression) {
		this.expression = expression;
	}

	public Gene(int geneId, List<Double> expression) {
		this.geneId = geneId;
		this.expression = expression;
	}
}
