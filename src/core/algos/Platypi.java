package core.algos;

import java.io.BufferedWriter;
import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.util.ArrayList;
import java.util.List;

import com.google.common.collect.EvictingQueue;

import core.algos.PlacementProcess;
import core.algos.SequenceKnife;
import inputs.Fasta;
import it.unimi.dsi.fastutil.chars.Char2FloatMap;
import main_v2.SessionNext_v2;

public class Platypi {
	class NodeStats {
		// tab[#times_encoutered] --> index=nodeId
		int[] nodeOccurences;
		// tab[score] --> index=nodeId
		float[] nodeScores;
		// copy that will be consumed in the Hoare's selection algorithm
		float[] nodeScoresCopy; 
		
		PlacementProcess.Score[] bestScoreList;
		
		public NodeStats(SessionNext_v2 session, int keepAtMost) {
			nodeOccurences = new int[session.originalTree.getNodeCount()];
			nodeScores = new float[session.originalTree.getNodeCount()];
			nodeScoresCopy = new float[session.originalTree.getNodeCount()]; 
			bestScoreList = new PlacementProcess.Score[keepAtMost];
		}
		
		public void init() {
			for (int i = 0; i < bestScoreList.length; i++) {
				bestScoreList[i] = new PlacementProcess.Score(-1, Float.NEGATIVE_INFINITY);
			}
		}
	}
	
	class History {
		int[] nodeOccurences;
		float[] nodeScores;
		
		public History(SessionNext_v2 session) {
			nodeOccurences = new int[session.originalTree.getNodeCount()];
			nodeScores = new float[session.originalTree.getNodeCount()];
		}
	}
	
	public static class Score {
		int nodeId;
		double weightRatio;
		
		public Score(int nodeId, double weightRatio) {
			this.nodeId = nodeId;
			this.weightRatio = weightRatio;
		}
	}
	
	public static class Placement {
		List<Score> global;
		List<List<Score>> windows;
		
		public Placement(List<Score> global, List<List<Score>> windows) {
			this.global = global;
			this.windows = windows;
		}
	}
	
	private NodeStats global;
	private NodeStats window;
	
	private SessionNext_v2 session;
	private int keepAtMost;
	
	public Platypi(SessionNext_v2 session, int keepAtMost) {
		this.session = session;
		this.keepAtMost = keepAtMost;
		
		global = new NodeStats(session, keepAtMost);
		window = new NodeStats(session, keepAtMost);
	}
	
	public Placement place(Fasta fasta, File logDir, int s, int w) throws IOException {
		boolean merStats = false;
		
		global.init();
		window.init();

		BufferedWriter bwMerStats = null;
		if (merStats) {
			File fileMerStats = new File(logDir + File.separator + "merStats.csv");
			bwMerStats = Files.newBufferedWriter(fileMerStats.toPath());
			PlacementProcess.merStatsHeader(bwMerStats);
		}

		SequenceKnife sk = new SequenceKnife(fasta, session.k, session.minK, session.states, SequenceKnife.SAMPLING_LINEAR);

		//TODO: room for optimization
		EvictingQueue<History> windowQueue = EvictingQueue.create(w);
		
		List<List<Score>> windowScores = new ArrayList<List<Score>>((fasta.getSequence(false).length()/w)+1);
		
		// merFound[nodeId][merPos]
		boolean[][] merFound = new boolean[session.ARTree.getNodeCount()][sk.getMerCount()];
		// loop on words
		byte[] qw = null;
		//TODO: queryKmerCount is just the amount of positions considered?
		int pos = 0;
		
		while ((qw = sk.getNextByteWord()) != null) {
			// get Pairs associated to this word
			Char2FloatMap.FastEntrySet allPairs = session.hash.getPairsOfTopPosition2(qw);

			History h = new History(session);
			
			// stream version, 5-10% faster than allPairs.fastIterator()
			allPairs.stream().forEach((Char2FloatMap.Entry entry) -> {
				int nodeId = entry.getCharKey();
				// count # times originalNode encountered
				global.nodeOccurences[nodeId] += 1;
				// score associated to originalNode x for current read
				global.nodeScores[nodeId] += entry.getFloatValue();
				
				window.nodeOccurences[nodeId] += 1;
				window.nodeScores[nodeId] += entry.getFloatValue();
				
				h.nodeOccurences[nodeId] += 1;
				h.nodeScores[nodeId] += entry.getFloatValue();				
			});
			
			windowQueue.add(h);
			
			//we reached a slice of the sliding window
			if (pos >= w && (pos + 1) % s == 0) {
				//score the slice
				List<Score> scores = score(window, keepAtMost);
				windowScores.add(scores);
				
				//shift the window
				for (int i = 0; i < session.originalTree.getNodeCount(); i++) {
					for (int j = 0; j < s; j++) {
						History hh = windowQueue.poll();
						window.nodeOccurences[i] -= hh.nodeOccurences[i];
						window.nodeScores[i] -= hh.nodeScores[i];
					}
				}
			}
			
			pos++;
		}
		
		if (merStats) {
			PlacementProcess.merStats(session, merFound, fasta, bwMerStats);
		}
		
		List<Score> scores = score(global, keepAtMost);

		return new Placement(scores, windowScores);
	}
	
	private static List<Score> score(NodeStats stats, int keepAtMost) {
		//TODO: take into account keepFactor
		
		//TODO: define selectedNodes to be all the nodes
		//TODO: we need to keep the selectedNodes explicitly
		ArrayList<Integer> selectedNodes_ = new ArrayList<Integer>();
		for (int n = 0; n < stats.nodeOccurences.length; n++) {
			selectedNodes_.add(n);
		}
		int numberOfBestScoreToConsiderForOutput 
			= PlacementProcess.selectionAlgo(
			                		keepAtMost,
			                		stats.nodeScores, stats.nodeScoresCopy, 
			                		stats.bestScoreList,
			                		selectedNodes_);
            
        double weightRatioShift = PlacementProcess.computeWeightRatioShift(stats.bestScoreList[0]);
        double allLikelihoodSums = PlacementProcess.computeLikelihoodSum(stats.bestScoreList, numberOfBestScoreToConsiderForOutput, weightRatioShift);
	
        List<Score> scores = new ArrayList<Score>();
        for (int i = stats.bestScoreList.length-1; i>stats.bestScoreList.length-numberOfBestScoreToConsiderForOutput-1; i--) {
            double weigth_ratio=PlacementProcess.computeWeightRatio(stats.bestScoreList[i], weightRatioShift, allLikelihoodSums);
            scores.add(new Score(stats.bestScoreList[i].nodeId, weigth_ratio));
        }
        return scores;
	}
}
