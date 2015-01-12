/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package dm_p2;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.Scanner;
import org.apache.commons.math3.stat.correlation.PearsonsCorrelation;

/**
 *
 * @author Deepak
 */
public class DBSCAN {

    public static class cluster {

        private int clid;
        private ArrayList<gene> genelist;

        cluster(int id) {
            clid = id;
            genelist = new ArrayList<gene>();
        }

        public void addtocl(gene g) {
            genelist.add(g);
        }

        public ArrayList<gene> getlist() {
            return genelist;
        }

        public int getid() {
            return clid;
        }
    }

    public static class gene {

        private int gid;
        private ArrayList<Double> exp;
        private String cl;
        private int visited;
        private int clusid;

        gene() {
            visited = 0;
            cl = "";
            clusid = 0;
            exp = new ArrayList<Double>();
        }

        public void setclusid(int cid) {
            clusid = cid;
        }

        public int getclusid() {
            return clusid;
        }

        public void add(int g, ArrayList<Double> e) {
            gid = g;
            exp.addAll(e);
        }

        public void addcl(String s) {
            cl = s;
        }

        public String getcl() {
            return cl;
        }

        public int getid() {
            return gid;
        }

        public int isVisited() {
            return visited;
        }

        public void Visited(int v) {
            visited = v;
        }

        public ArrayList<Double> getexp() {
            return exp;
        }
    }
    static LinkedHashMap<Integer, List<Double>> linkedHashMap = new LinkedHashMap<Integer, List<Double>>();
    static ArrayList<Double> expressionValues;
    public static Map<Integer, Integer> given = new HashMap<Integer, Integer>();
    public static Map<Integer, Integer> calculated = new HashMap<Integer, Integer>();
    static List<gene> genelist = new ArrayList<gene>();

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

    public static void normalize(LinkedHashMap<Integer, List<Double>> lhm) {
        Double[][] mydata = new Double[lhm.size()][lhm.get(1).size()];
        Double[] min = new Double[lhm.get(1).size()];
        Double[] max = new Double[lhm.get(1).size()];
        for (int i = 1; i <= lhm.size(); i++) {
            List<Double> eachrow = lhm.get(i);
            for (int k = 0; k < eachrow.size(); k++) {
                mydata[i - 1][k] = eachrow.get(k);
            }
        }
        for (int i = 0; i < lhm.get(1).size(); i++) {
            min[i] = 1000.0;
            max[i] = -1000.0;
            for (int k = 0; k < lhm.size(); k++) {
                if (mydata[k][i] < min[i]) {
                    min[i] = mydata[k][i];
                }
                if (mydata[k][i] > max[i]) {
                    max[i] = mydata[k][i];
                }
            }
        }
        for (int i = 1; i <= lhm.size(); i++) {
            for (int k = 0; k < lhm.get(1).size(); k++) {
                mydata[i - 1][k] = (mydata[i - 1][k] - min[k]) / (max[k] - min[k]);
            }
        }
        for (int i = 1; i <= lhm.size(); i++) {
            gene g = new gene();
            ArrayList<Double> expval = new ArrayList<Double>();
            for (int k = 0; k < lhm.get(1).size(); k++) {
                expval.add(mydata[i - 1][k]);
                //mydata[i-1][k]=(mydata[i-1][k]-min[k])/(max[k]-min[k]);
            }
            g.add(i, expval);
            linkedHashMap.put(i, expval);
            genelist.add(g);
        }

    }

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

    public static double[][] distmat(Map<Integer, List<Double>> lhm) {
        System.out.println(lhm.size());
        double[][] dmat = new double[lhm.size() + 1][lhm.size() + 1];
        for (int i = 1; i <= lhm.size(); i++) {
            List<Double> curg = lhm.get(i);
            for (int j = 1; j <= lhm.size(); j++) {
                if (i != j) {
                    List<Double> nxtg = lhm.get(j);
                    double sum = 0;
                    for (int k = 0; k < curg.size(); k++) {
                        sum += Math.pow((curg.get(k) - nxtg.get(k)), 2);
                    }
                    dmat[i][j] = Math.sqrt(sum);
                } else {
                    dmat[i][j] = 0.0;
                }
            }
        }
        return dmat;
    }

    public static double kdist(gene g1, gene g2) {
        ArrayList<Double> gene1 = g1.getexp();
        ArrayList<Double> gene2 = g2.getexp();
        double distance = 0.0;
        for (int k = 0; k < gene1.size(); k++) {
            distance += Math.pow((gene1.get(k) - gene2.get(k)), 2);
        }
        return Math.sqrt(distance);
    }

    public static ArrayList<cluster> DBSCAN(double eps, int MinPts) {
        ArrayList<cluster> clist = new ArrayList<cluster>();
        int c = 1;
        for (int i = 0; i < genelist.size(); i++) {
            if (genelist.get(i).isVisited() == 0) {
                genelist.get(i).Visited(1);
                ArrayList<gene> npts = regionQuery(genelist.get(i), eps);
                if (npts.size() < MinPts) {
                    genelist.get(i).addcl("NOISE");
                } else {
                    cluster cl = new cluster(c);
                    genelist.get(i).setclusid(c);
                    cl.addtocl(genelist.get(i));
                    for (int k = 0; k < npts.size(); k++) {
                        if (npts.get(k).isVisited() == 0 || (npts.get(k).isVisited() == 1 && npts.get(k).getcl().compareTo("NOISE") == 0)) {
                            npts.get(k).Visited(1);
                            ArrayList<gene> npts1 = regionQuery(npts.get(k), eps);
                            if (npts1.size() >= MinPts) {
                                npts.addAll(npts1);
                            }
                            if (npts.get(k).getclusid() == 0) {
                                npts.get(k).setclusid(c);
                                cl.addtocl(npts.get(k));
                            }
                        }
                    }
                    clist.add(cl);
                    c++;

                }
            }
        }
        return clist;
    }

    public static ArrayList<gene> regionQuery(gene g, double eps) {
        ArrayList<gene> neighpts = new ArrayList<gene>();
        for (int i = 0; i < genelist.size(); i++) {
            double dist2 = kdist(g, genelist.get(i));
            if (dist2 < eps) {
                neighpts.add(genelist.get(i));
            }
        }
        return neighpts;

    }

    public static void internalValidation(LinkedHashMap<Integer, List<Double>> linkedMap, Map<Integer, Integer> ourMap) {
        int size = linkedMap.size();
        int n = ((size - 1) * (size)) / 2;

        double incidenceMatrix[] = new double[n];
        double distanceMatrix[] = new double[n];
        int k = 0;
        for (int i = 1; i <= size; i++) {
            for (int j = i + 1; j <= size; j++) {
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

    public void run_dbscan() {
        generateLinkedHashMap("cho.txt");
        normalize(linkedHashMap);
        ArrayList<cluster> myclist = DBSCAN(0.257, 5);
        for (int i = 0; i < myclist.size(); i++) {
            cluster temp = myclist.get(i);
            ArrayList<gene> mycgene = temp.getlist();
            for (int k = 0; k < mycgene.size(); k++) {
                calculated.put(mycgene.get(k).getid(), temp.getid());
            }
        }
        externalValidation(given, calculated);
        internalValidation(linkedHashMap, calculated);

    }

}
