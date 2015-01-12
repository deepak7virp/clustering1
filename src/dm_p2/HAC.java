/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package dm_p2;
import static dm_p2.DBSCAN.expressionValues;
import static dm_p2.DBSCAN.generateLinkedHashMap;
import java.io.File;
import java.io.FileNotFoundException;
import java.util.*;
/**
 *
 * @author Deepak
 */
public class HAC {

    public Map<Integer, Integer> ourMap = new HashMap<Integer, Integer>();
    public Map<Integer, ArrayList<Integer>> clusters = new HashMap<Integer, ArrayList<Integer>>();
    public static Map<Integer, Integer> given = new HashMap<Integer, Integer>();
    public static Map<Integer, Integer> calculated = new HashMap<Integer, Integer>();
    public Map<Integer, Map<Integer, ArrayList<Integer>>> dendogram = new TreeMap<Integer, Map<Integer, ArrayList<Integer>>>();
    static LinkedHashMap<Integer, List<Double>> linkedHashMap = new LinkedHashMap<Integer, List<Double>>();
    public static void generateLinkedHashMap(String filename) {
        String filePath = new File("").getAbsolutePath();

        try {
            Scanner s = new Scanner(new File(filePath + "/src/dm_p2/" + filename));
            while (s.hasNext()) {
                String inputLine = s.nextLine();
                String[] splitData = inputLine.split("\t");
                int geneId = Integer.parseInt(splitData[0]);
                expressionValues = new ArrayList<Double>();
                int cc = Integer.parseInt(splitData[1]);
                given.put(geneId, cc);
                for (int i = 2; i < splitData.length; i++) {

                    expressionValues.add(Double.parseDouble(splitData[i]));
                }
//                gene g=new gene();
//                g.add(geneId, expressionValues);
//                genelist.add(g);
                linkedHashMap.put(geneId, expressionValues);
            }
        } catch (FileNotFoundException e) {
            e.printStackTrace();
        }

    }

    public void haclustering(int cut) {
        generateLinkedHashMap("cho.txt");
        int size = linkedHashMap.size();
        double distMatrix[][] = new double[size+1][size+1];
        double distance;
        double minDist = Double.MAX_VALUE;
        int mini = 0, minj = 0;

        for (int i = 1; i <= size; i++) {
            for (int j = 1; j <= size; j++) {
                distance = euclidianDistance(linkedHashMap.get(i), linkedHashMap.get(j));
                if (i == j) {
                    distMatrix[i][j] = Double.MAX_VALUE;
                } else {
                    distMatrix[i][j] = distMatrix[j][i] = distance;
                }

                if (distMatrix[i][j] < minDist) {
                    minDist = distMatrix[i][j];
                    mini = i;
                    minj = j;
                }
            }
        }

        //Initial Cluster to Gene (1:1) Mapping  
        Iterator itr = linkedHashMap.keySet().iterator();
        int clusterCount = 0;
        while (itr.hasNext() && clusterCount <= size) {
            Integer geneId = (Integer) itr.next();
            ArrayList<Integer> temp = new ArrayList<Integer>();
            temp.add(geneId);
            clusters.put(clusterCount, temp);
            clusterCount++;
        }

        int level = clusters.size();
        while (clusters.size() != 1) {
            if (mini < minj) {
                distMatrix = combine(mini, minj, clusters, distMatrix);
            } else if (minj < mini) {
                distMatrix = combine(minj, mini, clusters, distMatrix);
            } else {
                distMatrix = combine(mini, minj, clusters, distMatrix);
            }
            dendogram.put(level, clusters);
            level--;
        }

        for (int a = 0; a < clusters.size(); a++) {
            List<Integer> geneIdList = clusters.get(a);
            for (int b = 0; b < geneIdList.size(); b++) {
                int geneId = geneIdList.get(b);
                ourMap.put(geneId, a + 1);
            }
        }
    }

    public Map<Integer, ArrayList<Integer>> clusterInfo(int cutlevel) {
        Map<Integer, ArrayList<Integer>> clusterAtLevel = dendogram.get(cutlevel);
        return clusterAtLevel;
    }
    
    public double[][] combine(int small, int big, Map<Integer, ArrayList<Integer>> clusters, double distMatrix[][]) {
        int size = clusters.size();
        
        if (small < big) {
            clusters.get(small-1).addAll(clusters.get(big-1));
            clusters.remove(big-1);
            for (int i = 1; i <= size; i++) {
                if (distMatrix[i][big] < distMatrix[i][small]) {
                    distMatrix[i][small] = distMatrix[i][big];
                    distMatrix[small][i] = distMatrix[big][i];
                }
                distMatrix[i][big] = Double.MAX_VALUE;
                distMatrix[big][i] = Double.MAX_VALUE;
            }
        } else if (small == big) {
            distMatrix[small][small] = Double.MAX_VALUE;
        }
        return distMatrix;
    }

    public double euclidianDistance(List<Double> listi, List<Double> listj) {
        double euclidianDistance = 0;
        for (int i = 0; i < listi.size(); i++) {
            double distance = listi.get(i) - listj.get(i);
            euclidianDistance += Math.pow(distance, 2);
        }
        return Math.sqrt(euclidianDistance);
    }
}
