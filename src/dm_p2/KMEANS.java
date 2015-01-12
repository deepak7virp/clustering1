/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package dm_p2;

import static dm_p2.DBSCAN.linkedHashMap;
import java.io.File;
import java.io.FileNotFoundException;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.Scanner;
import java.util.TreeMap;
import org.apache.commons.math3.stat.correlation.PearsonsCorrelation;

/**
 *
 * @author Deepak
 */
public class KMEANS {

    public static class Gene {

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
    static Map<Integer, List<Double>> linkedHashMap = new LinkedHashMap<Integer, List<Double>>();
    static List<Double> expressionValues;
    static List<Gene> centroids = new ArrayList<>();
    static List<List<Gene>> clusters = new ArrayList<List<Gene>>();
    static int[] initialCentroidsRows = {1, 7, 9, 10, 15};
    static TreeMap<Integer, Integer> geneToCluster = new TreeMap<>();
    static TreeMap<Integer, Integer> validationMap = new TreeMap<>();
    public static void externalValidation(Map<Integer, Integer> givenMap, Map<Integer, Integer> ourMap) {
        int truthMatrix[][] = new int[givenMap.size()][givenMap.size()];
        int clusterMatrix[][] = new int[givenMap.size()][givenMap.size()];

        int SS = 0, SD = 0, DS = 0, DD = 0;
        for (int i = 0; i < truthMatrix.length; i++) {
            for (int j = 0; j < clusterMatrix.length; j++) {
                if (givenMap.get(i) == givenMap.get(j)) {
                    truthMatrix[i][j] = truthMatrix[j][i] = 1;
                } else {
                    truthMatrix[i][j] = truthMatrix[j][i] = 0;
                }
                if (ourMap.get(i) == ourMap.get(j)) {
                    clusterMatrix[i][j] = clusterMatrix[j][i] = 1;
                } else {
                    clusterMatrix[i][j] = clusterMatrix[j][i] = 0;
                }
                if (truthMatrix[i][j] == clusterMatrix[i][j] && truthMatrix[i][j] == 1) {
                    SS++;
                } else if (truthMatrix[i][j] == clusterMatrix[i][j] && truthMatrix[i][j] == 0) {
                    DD++;
                } else if ((truthMatrix[i][j] != clusterMatrix[i][j]) && (truthMatrix[i][j] == 0 && clusterMatrix[i][j] == 1)) {
                    DS++;
                } else if ((truthMatrix[i][j] != clusterMatrix[i][j]) && (truthMatrix[i][j] == 1 && clusterMatrix[i][j] == 0)) {
                    SD++;
                }
            }
        }
        System.out.println("External Index Validation : Rand Index = " + (SS + DD) / (double) (SS + SD + DS + DD));
        System.out.println("External Index Validation : Jaccard Coefficient = " + (SS) / (double) (SS + SD + DS));
    }
    public static void generateLinkedHashMap(String filename) {
        String filePath = new File("").getAbsolutePath();
        try {
            Scanner s = new Scanner(new File(filePath + "/src/dm_p2/" + filename));
            while (s.hasNext()) {
                String inputLine = s.nextLine();
                String[] splitData = inputLine.split("\t");
                int geneId = Integer.parseInt(splitData[0]);
                int valgeneId = Integer.parseInt(splitData[1]);
                expressionValues = new ArrayList<Double>();
                for (int i = 2; i < splitData.length; i++) {
                    expressionValues.add(Double.parseDouble(splitData[i]));
                }
                linkedHashMap.put(geneId, expressionValues);
                if (valgeneId != -1) {
                    validationMap.put(geneId, valgeneId);
                }
            }
            //System.out.println(linkedHashMap.entrySet());
            //System.out.println(validationMap.entrySet());
        } catch (FileNotFoundException e) {
            e.printStackTrace();
        }
    }
    public static void internalValidation(Map<Integer, List<Double>> linkedMap, Map<Integer, Integer> ourMap) {
        int size = linkedMap.size();
        int n = ((size-1) * (size)) / 2;

        double incidenceMatrix[] = new double[n];
        double distanceMatrix[] = new double[n];
        int k=0;
        for (int i = 1; i <= size; i++) {
            for (int j = i+1; j <= size; j++) {
                if (ourMap.get(i) == ourMap.get(j)) {
                    incidenceMatrix[k] = 1;
                } else {
                    incidenceMatrix[k] = 0;
                }
                
                distanceMatrix[k] = euclidianDistance(linkedMap.get(i), linkedMap.get(j));
                k++;
            }
        }
        PearsonsCorrelation pc = new PearsonsCorrelation();
        
        System.out.println("Internal Index Validation : Correlation = " + pc.correlation(incidenceMatrix, distanceMatrix));

    }

    public static double euclidianDistance(List<Double> listi, List<Double> listj) {
        double euclidianDistance = 0;
        for (int i = 0; i < listi.size(); i++) {
            double distance = listi.get(i) - listj.get(i);
            euclidianDistance += Math.pow(distance, 2);
        }
        return Math.sqrt(euclidianDistance);
    }
    
	// First defining a utility functions needed for K Means
    //Calculating mean values (Centroids) of clusters
    //Input is a list of gene and each gene has a geneid and list of expression values associated to it
    public static ArrayList<Double> calculateCentroids(List<Gene> input) {
        ArrayList<Double> avgList = new ArrayList<>();
        double avg, sum;
        for (int i = 0; i < input.get(0).getExpression().size(); i++) {
            sum = 0;
            avg = 0;
            for (int j = 0; j < input.size(); j++) {
                sum += input.get(j).getExpression().get(i);
            }
            avg = sum / input.size();
            avgList.add(avg);
        }
        return avgList;
    }

    public static void normalize(Map<Integer, List<Double>> lhm){
        Double[][] mydata=new Double[lhm.size()][lhm.get(1).size()];
        Double[] min=new Double[lhm.get(1).size()];
        Double[] max=new Double[lhm.get(1).size()];
        for(int i=1;i<=lhm.size();i++){
            List<Double> eachrow=lhm.get(i);
            for(int k=0;k<eachrow.size();k++){
                mydata[i-1][k]=eachrow.get(k);
            }
        }
        for(int i=0;i<lhm.get(1).size();i++){
            min[i]=1000.0;
            max[i]=-1000.0;
            for(int k=0;k<lhm.size();k++){
                if(mydata[k][i]<min[i]){
                    min[i]=mydata[k][i];
                }
                if(mydata[k][i]>max[i]){
                    max[i]=mydata[k][i];
                }
            }
        }
        for(int i=1;i<=lhm.size();i++){
            for(int k=0;k<lhm.get(1).size();k++){
                mydata[i-1][k]=(mydata[i-1][k]-min[k])/(max[k]-min[k]);
            }
        }
        for(int i=1;i<=lhm.size();i++){
            //DBSCAN.gene g=new DBSCAN.gene();
            ArrayList<Double> expval=new ArrayList<Double>();
            for(int k=0;k<lhm.get(1).size();k++){
                expval.add(mydata[i-1][k]);
                //mydata[i-1][k]=(mydata[i-1][k]-min[k])/(max[k]-min[k]);
            }
            //g.add(i, expval);
            linkedHashMap.put(i, expval);
            //genelist.add(g);
        }
        
        
    }
    
    //K-Means Algorithm
    public static void runKMeans(Map<Integer, List<Double>> linkedHashMap, int[] initialCentroidsRows, int numOfClusters, int numOfIterations) {

        int clusterCount = 0;
//        List<Gene> geneData;
//        for (int i = 0; i < numOfClusters; i++) {
//            geneData = new ArrayList<Gene>();
//            clusters.add(geneData);
//        }
        for (int k = 0; k < numOfClusters; k++) {
            Integer geneId = initialCentroidsRows[clusterCount++];
            Gene gene = new Gene(geneId, linkedHashMap.get(geneId));
            centroids.add(gene);
            //clusterCount++;
        }
        
        for (int num = 0; num < numOfIterations; num++) {
            Iterator itr = linkedHashMap.keySet().iterator();
            clusterCount=0;
            List<Gene> geneData1;
            for (int i = 0; i < numOfClusters; i++) {
                geneData1 = new ArrayList<Gene>();
                clusters.add(geneData1);
            }
            while (itr.hasNext()) {
                double minDist = Double.MAX_VALUE;
                int minDistIndex = 0;
                Integer geneId = (Integer) itr.next();

                Gene gene = new Gene(geneId, linkedHashMap.get(geneId));
                List<Double> geneExpressionValues = new ArrayList<Double>();
                geneExpressionValues = linkedHashMap.get(geneId);
                for (int j = 0; j < centroids.size(); j++) {
                    double eucDistance = 0.0;
                    List<Double> centroidExpressionValues = new ArrayList<Double>();
                    centroidExpressionValues = centroids.get(j).getExpression();
                    for (int i = 0; i < centroidExpressionValues.size(); i++) {
                        double dist = geneExpressionValues.get(i) - centroidExpressionValues.get(i);
                        eucDistance += Math.pow(dist, 2);
                    }
                    double dist = Math.sqrt(eucDistance);
                    if (dist < minDist) {
                        minDist = dist;
                        minDistIndex = j;
                    }
                }
                //System.out.println(minDist+"-"+gene.geneId+"-"+centroids.size());
                clusters.get(minDistIndex).add(gene);

				// Storing the result of each iteration in a TreeMap needed for validation.
                // This TreeMap records the cluster to which each gene gets allocated 
                geneToCluster.put(geneId, minDistIndex);
            }
            centroids.clear();
            
             
            for (int i = 0; i < numOfClusters; i++) {
                ArrayList<Double> expValuesNew = new ArrayList<Double>();
                expValuesNew = calculateCentroids(clusters.get(i));
                Gene gene = new Gene(i, expValuesNew);
                centroids.add(i, gene);
            }
            if(num<numOfIterations-1)
            clusters.clear();
            
        }
    }

    public  void run() {
        generateLinkedHashMap("cho.txt");
        normalize(linkedHashMap);
        runKMeans(linkedHashMap, initialCentroidsRows, 5, 15);
        //System.out.println(geneToCluster.entrySet());
        internalValidation(linkedHashMap,geneToCluster);
        externalValidation(validationMap,geneToCluster);
    }

}
