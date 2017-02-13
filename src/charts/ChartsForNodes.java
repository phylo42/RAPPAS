/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package charts;

import alignement.Alignment;
import etc.Infos;
import grad.GradientSegment;
import java.awt.BasicStroke;
import java.awt.Color;
import java.awt.Dimension;
import java.awt.GradientPaint;
import java.awt.Graphics2D;
import java.awt.GridLayout;
import java.awt.Paint;
import java.awt.geom.Rectangle2D;
import java.util.Arrays;
import javax.swing.JFrame;
import javax.swing.JPanel;
import javax.swing.SwingUtilities;
import org.jfree.chart.ChartFactory;
import org.jfree.chart.ChartPanel;
import org.jfree.chart.JFreeChart;
import org.jfree.chart.axis.AxisLocation;
import org.jfree.chart.axis.NumberAxis;
import org.jfree.chart.axis.ValueAxis;
import org.jfree.chart.entity.EntityCollection;
import org.jfree.chart.plot.CategoryPlot;
import org.jfree.chart.plot.CrosshairState;
import org.jfree.chart.plot.PlotOrientation;
import org.jfree.chart.plot.PlotRenderingInfo;
import org.jfree.chart.plot.XYPlot;
import org.jfree.chart.renderer.PaintScale;
import org.jfree.chart.renderer.xy.XYBlockRenderer;
import org.jfree.chart.renderer.xy.XYItemRendererState;
import org.jfree.chart.title.PaintScaleLegend;
import org.jfree.data.category.DefaultCategoryDataset;
import org.jfree.data.xy.MatrixSeriesCollection;
import org.jfree.data.xy.XYDataset;
import org.jfree.data.xy.XYSeriesCollection;
import org.jfree.data.xy.XYZDataset;
import org.jfree.ui.RectangleEdge;
import org.jfree.ui.RectangleInsets;
import org.jfree.ui.RefineryUtilities;
import tree.PhyloNode;

/**
 *
 * @author ben
 */
public class ChartsForNodes extends javax.swing.JFrame {
    
private static final int N = 1000;


    /**
     * Creates new form NewJFrame
     */
    public ChartsForNodes() {
    }

    
    public static JPanel buildReadMatchForANode3(PhyloNode n, Alignment align, int readLength,XYZDataset dataset,double scaleLowerBound, double scaleUpperBound) {
        
        NumberAxis xAxis = new NumberAxis("Reference");
        xAxis.setUpperBound(align.getLength());
        xAxis.setLowerBound(0);
        NumberAxis yAxis = new NumberAxis("Read");
        yAxis.setUpperBound(readLength);
        yAxis.setLowerBound(0);
        
//        for (int i = 0; i < dataset.getItemCount(0); i++) {
//            System.out.println(dataset.getX(0, i)+","+dataset.getY(0, i)+","+dataset.getZ(0, i));
//        }
        
        
        
        XYPlot plot = new XYPlot(dataset, xAxis, yAxis, null);
        
//        XYBlockRenderer r = new XYBlockRenderer();
//        SpectrumPaintScale ps = new SpectrumPaintScale(0, 1);
//        r.setPaintScale(ps);
//        r.setBlockHeight(10.0f);
//        r.setBlockWidth(10.0f);
//        plot.setRenderer(r);
        
        //gray scale
        // OR full color scale !
        //PaintScale ps= new GrayPaintScale(scaleLowerBound, scaleUpperBound);
        //PaintScale ps= new SpectrumPaintScale(scaleLowerBound, scaleUpperBound);
        PaintScale ps= new MyPaintScale(scaleLowerBound, scaleUpperBound);
        
        
        
        XYBlockRenderer r = new XYBlockRenderer() {
            
            PaintScale paintScale= ps;
            double xOffset;
            double yOffset;
            double blockWidth = 1.0;
            double blockHeight = 1.0;
            
            @Override
            public void drawItem(Graphics2D g2, XYItemRendererState state, Rectangle2D dataArea, PlotRenderingInfo info, XYPlot plot, ValueAxis domainAxis, ValueAxis rangeAxis, XYDataset dataset, int series, int item, CrosshairState crosshairState, int pass) {
                double x = dataset.getXValue(series, item);
                double y = dataset.getYValue(series, item);
                double z = 0.0;
                if (dataset instanceof XYZDataset) {
                    z = ((XYZDataset) dataset).getZValue(series, item);
                }
                Paint p = this.paintScale.getPaint(z);
                double xx0 = domainAxis.valueToJava2D(x + this.xOffset, dataArea,
                        plot.getDomainAxisEdge());
                double yy0 = rangeAxis.valueToJava2D(y + this.yOffset, dataArea,
                        plot.getRangeAxisEdge());
                double xx1 = domainAxis.valueToJava2D(x + this.blockWidth
                        + this.xOffset, dataArea, plot.getDomainAxisEdge());
                double yy1 = rangeAxis.valueToJava2D(y + this.blockHeight
                        + this.yOffset, dataArea, plot.getRangeAxisEdge());
                Rectangle2D block;
                PlotOrientation orientation = plot.getOrientation();
                
                //modified below to double dot sizes
                if (orientation.equals(PlotOrientation.HORIZONTAL)) {
                    block = new Rectangle2D.Double(Math.min(yy0, yy1),
                            Math.min(xx0, xx1), 2*Math.abs(yy1 - yy0),
                            2*Math.abs(xx0 - xx1));
                }
                else {
                    block = new Rectangle2D.Double(Math.min(xx0, xx1),
                            Math.min(yy0, yy1), 2*Math.abs(xx1 - xx0),
                            2*Math.abs(yy1 - yy0));
                }
                //overidden to add transparency
                Color c= (Color) p;
                Color c1 = new Color(c.getRed(),c.getGreen(),c.getBlue(), 0); 
                g2.setPaint(p);
                //overidden to add transparency
                g2.fill(block);

                //g2.setStroke(new BasicStroke(5.0f));                    
                if (dataset.getY(series, item).doubleValue() > 0) {
                    g2.draw(block);
                }

                EntityCollection entities = state.getEntityCollection();
                if (entities != null) {
                    addEntity(entities, block, dataset, series, item, 0.0, 0.0);
                }            
            }
            
        };
        
        //r.setPaintScale(grayScale);   
//        r.setBlockHeight(1);
//        r.setBlockWidth(1);
        plot.setRenderer(r);
        
        JFreeChart chart = new JFreeChart("Mapping for nodeId="+n.getId()+" nodeName="+n.getLabel(),JFreeChart.DEFAULT_TITLE_FONT, plot, false);
        NumberAxis scaleAxis = new NumberAxis("log10(PP*)");
        scaleAxis.setAxisLinePaint(Color.white);
        scaleAxis.setTickMarkPaint(Color.white);
        
        //PaintScaleLegend legend = new PaintScaleLegend(ps, scaleAxis);
        PaintScaleLegend legend = new PaintScaleLegend(ps, scaleAxis);
        
        legend.setSubdivisionCount(128);
        legend.setAxisLocation(AxisLocation.TOP_OR_RIGHT);
        legend.setPadding(new RectangleInsets(10, 10, 10, 10));
        legend.setStripWidth(20);
        legend.setPosition(RectangleEdge.RIGHT);
        legend.setBackgroundPaint(Color.WHITE);
        chart.addSubtitle(legend);
        chart.setBackgroundPaint(Color.white);
        ChartPanel panel=new ChartPanel(chart);
        panel.setSize(1024, 250);
        panel.setPreferredSize(new Dimension(1024, 250));
        return panel;
        
    }
    
    public static JPanel buildReadMatchForANode2(PhyloNode n, Alignment align, int readLength,XYSeriesCollection dataset) {
        
        JFreeChart chart = ChartFactory.createScatterPlot(
                    "Mapping for nodeId="+n.getId()+" nodeName="+n.getLabel(), "Reference", "Read", dataset,
                    PlotOrientation.VERTICAL, true, true, false);

        chart.setBackgroundPaint(new GradientPaint(0, 0, Color.white, 0,1000, Color.blue));

        final XYPlot plot = chart.getXYPlot();
        plot.setForegroundAlpha(0.5f);

        final NumberAxis domainAxis = (NumberAxis) plot.getDomainAxis();
        domainAxis.setLowerBound(0);
        domainAxis.setUpperBound(align.getLength());
        final NumberAxis rangeAxis = (NumberAxis) plot.getRangeAxis();
        // rangeAxis.setInverted(true);  // uncoment to reproduce a bug in jFreeChart
        rangeAxis.setLowerBound(0);
        rangeAxis.setUpperBound(readLength);
        

        final ChartPanel chartPanel = new ChartPanel(chart);
        
        
        return chartPanel;
        
    }
    
    public static JPanel buildReadMatchForANode(PhyloNode n, MatrixSeriesCollection dataset) {
        
        final JFreeChart chart = ChartFactory.createBubbleChart(
            "TITLE", "X", "Y", dataset, 
            PlotOrientation.VERTICAL, 
            true,
            true, false);

        chart.setBackgroundPaint(new GradientPaint(0, 0, Color.white, 0,1000, Color.blue));

        final XYPlot plot = chart.getXYPlot();
        plot.setForegroundAlpha(0.5f);

        final NumberAxis domainAxis = (NumberAxis) plot.getDomainAxis();
        domainAxis.setLowerBound(0);
        final NumberAxis rangeAxis = (NumberAxis) plot.getRangeAxis();
        // rangeAxis.setInverted(true);  // uncoment to reproduce a bug in jFreeChart
        rangeAxis.setLowerBound(0);

        final ChartPanel chartPanel = new ChartPanel(chart);
        
        
        return chartPanel;
        
    }
     
    public static JPanel buildPPHistogramForANode(PhyloNode n, double[] ppValuesOfANode) {
        
        //for now, we will do 10 bars, for powers of 10 from 0 to 9
        
        int[] binSize=new int[10];
        Arrays.fill(binSize, 0, 9, 0);
        String[] categorie=new String[10];
        

        
        for (int i = 0; i < ppValuesOfANode.length; i++) {
            int scale = new Double(Math.abs(Math.floor(Math.log10(ppValuesOfANode[i])))).intValue();
            if (scale<10) {
                binSize[scale]=binSize[scale]+1;
            } else {
                binSize[9]=binSize[9]+1;
            }
        }
        int max=0;
        for (int i = 0; i < 9; i++) {
            categorie[i]="-"+Integer.toString(i);
            if (binSize[i]>max) {
                max=binSize[i];
            }
        }
        categorie[9]=">=-9";
        
        DefaultCategoryDataset dataset = new DefaultCategoryDataset();
        
        for (int i = 0; i < 10; i++) {
            dataset.addValue(binSize[i],"Bins", categorie[i]);
        }
        Infos.println("init JFreeChart: "+ppValuesOfANode.length+" probas");
        JFreeChart barChart = ChartFactory.createBarChart("pp values for node:\n"+n, "pp value (log10 bins) ", "# words in this node", dataset);
        
//        // OPTIONAL CUSTOMISATION 
//        // set the background color for the chart...
//        barChart.setBackgroundPaint(Color.white);
//
        // get a reference to the plot for further customisation...
        CategoryPlot plot = barChart.getCategoryPlot();
        plot.setBackgroundPaint(Color.lightGray);
        plot.setDomainGridlinePaint(Color.white);
        plot.setRangeGridlinePaint(Color.white);

        // set the range axis to display integers only...
        final NumberAxis rangeAxis = (NumberAxis) plot.getRangeAxis();
        rangeAxis.setStandardTickUnits(NumberAxis.createIntegerTickUnits());
        rangeAxis.setLowerBound(0);
        rangeAxis.setUpperBound(max);
//
//        // disable bar outlines...
//        final BarRenderer renderer = (BarRenderer) plot.getRenderer();
//        renderer.setDrawBarOutline(false);
//        
//        // set up gradient paints for series...
//        final GradientPaint gp0 = new GradientPaint(
//            0.0f, 0.0f, Color.blue, 
//            0.0f, 0.0f, Color.lightGray
//        );
//        final GradientPaint gp1 = new GradientPaint(
//            0.0f, 0.0f, Color.green, 
//            0.0f, 0.0f, Color.lightGray
//        );
//        final GradientPaint gp2 = new GradientPaint(
//            0.0f, 0.0f, Color.red, 
//            0.0f, 0.0f, Color.lightGray
//        );
//        renderer.setSeriesPaint(0, gp0);
//        renderer.setSeriesPaint(1, gp1);
//        renderer.setSeriesPaint(2, gp2);
//
//        final CategoryAxis domainAxis = plot.getDomainAxis();
//        domainAxis.setCategoryLabelPositions(
//            CategoryLabelPositions.createUpRotationLabelPositions(Math.PI / 6.0)
//        );
        // OPTIONAL CUSTOMISATION COMPLETED.
        
        
        ChartPanel chartPanel = new ChartPanel(barChart);
        chartPanel.setPreferredSize(new Dimension(500, 500));
        chartPanel.setVisible(true);
        return chartPanel;
}
    
    
    /**
     * This method is called from within the constructor to initialize the form.
     * WARNING: Do NOT modify this code. The content of this method is always
     * regenerated by the Form Editor.
     */
    @SuppressWarnings("unchecked")
    // <editor-fold defaultstate="collapsed" desc="Generated Code">//GEN-BEGIN:initComponents
    private void initComponents() {

        setDefaultCloseOperation(javax.swing.WindowConstants.EXIT_ON_CLOSE);

        javax.swing.GroupLayout layout = new javax.swing.GroupLayout(getContentPane());
        getContentPane().setLayout(layout);
        layout.setHorizontalGroup(
            layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGap(0, 400, Short.MAX_VALUE)
        );
        layout.setVerticalGroup(
            layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGap(0, 300, Short.MAX_VALUE)
        );

        pack();
    }// </editor-fold>//GEN-END:initComponents

    /**
     * @param args the command line arguments
     */
    public static void main(String args[]) {
        /* Set the Nimbus look and feel */
        //<editor-fold defaultstate="collapsed" desc=" Look and feel setting code (optional) ">
        /* If Nimbus (introduced in Java SE 6) is not available, stay with the default look and feel.
         * For details see http://download.oracle.com/javase/tutorial/uiswing/lookandfeel/plaf.html 
         */
        try {
            for (javax.swing.UIManager.LookAndFeelInfo info : javax.swing.UIManager.getInstalledLookAndFeels()) {
                if ("Nimbus".equals(info.getName())) {
                    javax.swing.UIManager.setLookAndFeel(info.getClassName());
                    break;
                }
            }
        } catch (ClassNotFoundException ex) {
            java.util.logging.Logger.getLogger(ChartsForNodes.class.getName()).log(java.util.logging.Level.SEVERE, null, ex);
        } catch (InstantiationException ex) {
            java.util.logging.Logger.getLogger(ChartsForNodes.class.getName()).log(java.util.logging.Level.SEVERE, null, ex);
        } catch (IllegalAccessException ex) {
            java.util.logging.Logger.getLogger(ChartsForNodes.class.getName()).log(java.util.logging.Level.SEVERE, null, ex);
        } catch (javax.swing.UnsupportedLookAndFeelException ex) {
            java.util.logging.Logger.getLogger(ChartsForNodes.class.getName()).log(java.util.logging.Level.SEVERE, null, ex);
        }
        //</editor-fold>
        //</editor-fold>
        //</editor-fold>
        //</editor-fold>
        //</editor-fold>
        //</editor-fold>
        //</editor-fold>
        //</editor-fold>

        /* Create and display the form */
        SwingUtilities.invokeLater(new Runnable() {
            @Override
            public void run() {
                
                double[] test = {
                               1e-6,
                               1.24e-5,
                               1.235e-2,
                               1.95e-2,
                               0.245,
                               1.399e-11
                };
                
                
                ChartsForNodes demo=new ChartsForNodes();
                demo.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
                demo.setLayout(new GridLayout(1, 1));
                demo.add(demo.buildPPHistogramForANode(new PhyloNode(0, "TestNode:N0",0.1f), test));
                demo.pack();
                RefineryUtilities.centerFrameOnScreen(demo);
                demo.setVisible(true);
            }
        });
    }

    // Variables declaration - do not modify//GEN-BEGIN:variables
    // End of variables declaration//GEN-END:variables

    private static class SpectrumPaintScale implements PaintScale {

        private static final float H1 = 0f;
        private static final float H2 = 0.5f;  //SET HERE to 0.1f to get full scope
        private final double lowerBound;
        private final double upperBound;

        public SpectrumPaintScale(double lowerBound, double upperBound) {
            this.lowerBound = lowerBound;
            this.upperBound = upperBound;
        }

        @Override
        public double getLowerBound() {
            return lowerBound;
        }

        @Override
        public double getUpperBound() {
            return upperBound;
        }

        @Override
        public Paint getPaint(double value) {
            float scaledValue = (float) (value / (getUpperBound() - getLowerBound()));
            float scaledH = H1 + scaledValue * (H2 - H1);
            return Color.getHSBColor(scaledH, 1f, 0.9f);
        }
    }

    private static class MyPaintScale implements PaintScale {

           private final double lowerBound;
           private final double upperBound;

           public MyPaintScale(double lowerBound, double upperBound) {
               this.lowerBound = lowerBound;
               this.upperBound = upperBound;
           }

           @Override
           public double getLowerBound() {
               return lowerBound;
           }

           @Override
           public double getUpperBound() {
               return upperBound;
           }

           @Override
           public Paint getPaint(double value) {
               GradientSegment gs=new GradientSegment(new Color(200, 0, 0, 50), new Color(0, 0, 200, 50), 0.0,1.0,false);
               //scales on (-1;1)
               float scaledValue = (float) (value / (getUpperBound() - getLowerBound())); 
               if (scaledValue<=-1.0) {
                   return new Color(255, 255, 255, 255); //fully transparent, so invisible
               } else {
                   return gs.getColor(Math.abs(scaledValue));
               }
           }
    }
}
