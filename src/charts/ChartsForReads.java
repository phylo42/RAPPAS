/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package charts;

import java.awt.Color;
import java.awt.Dimension;
import java.awt.Graphics2D;
import java.awt.Paint;
import java.awt.geom.Rectangle2D;
import javax.swing.JPanel;
import org.jfree.chart.ChartFactory;
import org.jfree.chart.ChartPanel;
import org.jfree.chart.JFreeChart;
import org.jfree.chart.axis.NumberAxis;
import org.jfree.chart.axis.ValueAxis;
import org.jfree.chart.entity.EntityCollection;
import org.jfree.chart.plot.CrosshairState;
import org.jfree.chart.plot.PlotOrientation;
import org.jfree.chart.plot.PlotRenderingInfo;
import org.jfree.chart.plot.XYPlot;
import org.jfree.chart.renderer.PaintScale;
import org.jfree.chart.renderer.xy.XYBlockRenderer;
import org.jfree.chart.renderer.xy.XYItemRendererState;
import org.jfree.chart.renderer.xy.XYLineAndShapeRenderer;
import org.jfree.data.xy.XYDataset;
import org.jfree.data.xy.XYZDataset;
import tree.PhyloTree;

/**
 *
 * @author ben
 */
public class ChartsForReads {
    
    public static JPanel buildDiagScoreGraph(PhyloTree tree, XYDataset dataset, String read) {
        JFreeChart chart = ChartFactory.createXYLineChart(
            "Placement score for "+read,      // chart title
            "Reference",                      // x axis label
            "Score",                      // y axis label
            dataset,                  // data
            PlotOrientation.VERTICAL,
            true,                     // include legend
            true,                     // tooltips
            false                     // urls
        );

        // NOW DO SOME OPTIONAL CUSTOMISATION OF THE CHART...
        chart.setBackgroundPaint(Color.white);
        
        
        XYPlot plot = chart.getXYPlot();
        
        plot.setBackgroundPaint(Color.lightGray);
        plot.setDomainGridlinePaint(Color.white);
        plot.setRangeGridlinePaint(Color.white);
        
        //set the start/end of the X axis 100 bp above the reference
        plot.getDomainAxis().setLowerBound(-100);
        plot.getDomainAxis().setUpperBound(dataset.getItemCount(0)+100);
        
        
        XYLineAndShapeRenderer renderer = new XYLineAndShapeRenderer() {
            double xOffset;
            double yOffset;
            double blockWidth = 1.0;
            double blockHeight = 1.0;
            private boolean drawSeriesLineAsPath;

            @Override
            public void drawItem(Graphics2D g2, XYItemRendererState state, Rectangle2D dataArea, PlotRenderingInfo info, XYPlot plot, ValueAxis domainAxis, ValueAxis rangeAxis, XYDataset dataset, int series, int item, CrosshairState crosshairState, int pass) {
                // do nothing if item is not visible
                if (!getItemVisible(series, item)) {
                    return;
                }
                // first pass draws the background (lines, for instance)
                if (isLinePass(pass)) {
                    if (getItemLineVisible(series, item)) {
                        if (this.drawSeriesLineAsPath) {
                            //if (dataset.getY(series, item).doubleValue()!=0.0)
                                drawPrimaryLineAsPath(state, g2, plot, dataset, pass,
                                    series, item, domainAxis, rangeAxis, dataArea);
                        }
                        else {
                            //if (dataset.getY(series, item).doubleValue()!=0.0)
                                drawPrimaryLine(state, g2, plot, dataset, pass, series,
                                    item, domainAxis, rangeAxis, dataArea);
                        }
                    }
                }
                // second pass adds shapes where the items are ..
                else if (isItemPass(pass)) {
                    // setup for collecting optional entity info...
                    EntityCollection entities = null;
                    if (info != null && info.getOwner() != null) {
                        entities = info.getOwner().getEntityCollection();
                    }
                    if (dataset.getY(series, item).doubleValue()!=0.0)
                        drawSecondaryPass(g2, plot, dataset, pass, series, item,
                            domainAxis, dataArea, rangeAxis, crosshairState, entities);
                }
            }
        };
        renderer.setDrawSeriesLineAsPath(true);
        for (int i=0;i<plot.getSeriesCount();i++) {
            renderer.setSeriesShapesVisible(i, true);
            renderer.setSeriesLinesVisible(i, false);
        }
        plot.setRenderer(renderer);

        // change the auto tick unit selection to integer units only...
        final NumberAxis rangeAxis = (NumberAxis) plot.getRangeAxis();
        rangeAxis.setStandardTickUnits(NumberAxis.createIntegerTickUnits());
        // OPTIONAL CUSTOMISATION COMPLETED.
                
        ChartPanel panel=new ChartPanel(chart);
        panel.setSize(1024, 250);
        panel.setPreferredSize(new Dimension(1024, 250));
        return panel;
    }
    
}
