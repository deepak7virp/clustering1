/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package dm_p2;

/**
 *
 * @author Deepak
 */
public class DM_P2 {

    /**
     * @param args the command line arguments
     */
    public static void main(String[] args) {
        // TODO code application logic here
        KMEANS km = new KMEANS();
        System.out.println("__________K-MEANS_________");
        km.run();
        DBSCAN dbs = new DBSCAN();
        System.out.println("__________DBSCAN_________");
        dbs.run_dbscan();
        HAC hac=new HAC();
        hac.haclustering(0);
    }

}
