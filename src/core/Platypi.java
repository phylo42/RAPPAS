package core;

import java.io.BufferedWriter;
import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.util.ArrayList;

import core.algos.PlacementProcess;
import core.algos.SequenceKnife;
import inputs.Fasta;
import it.unimi.dsi.fastutil.chars.Char2FloatMap;
import main_v2.SessionNext_v2;

public class Platypi {
	public static void place(SessionNext_v2 session, Fasta fasta, int queryWordSampling, int keepAtMost, File logDir) {
		boolean merStats = false;

		ArrayList<Integer> selectedNodes = new ArrayList<>(10);
		// tab[#times_encoutered] --> index=nodeId
		int[] nodeOccurences = new int[session.originalTree.getNodeCount()];
		// tab[score] --> index=nodeId
		float[] nodeScores = new float[session.originalTree.getNodeCount()];
		// copy that will be consumed in the Hoare's selection algorithm
		float[] nodeScoresCopy = new float[session.originalTree.getNodeCount()];

		// System.out.println("S/C size: "+nodeOccurences.length);
		float kthLargestValue = -1;
		int numberOfBestScoreToConsiderForOutput = -1;
		Score[] bestScoreList = new Score[keepAtMost]; // in ascending order
		for (int i = 0; i < bestScoreList.length; i++) {
			bestScoreList[i] = new Score(-1, Float.NEGATIVE_INFINITY);
		}

		BufferedWriter bwMerStats = null;
		if (merStats) {
			File fileMerStats = new File(logDir + File.separator + "merStats.csv");
			try {
				bwMerStats = Files.newBufferedWriter(fileMerStats.toPath());
				PlacementProcess.merStatsHeader(bwMerStats);
			} catch (IOException e) {
				e.printStackTrace();
			}
		}

		SequenceKnife sk = new SequenceKnife(fasta, session.k, session.minK, session.states, queryWordSampling);
		int queryKmerCount = 0;
		int queryKmerMatchingDB = 0;

		// TODO: perhaps this loop can be shared between vanilla rappas and
		// platypi
		// merFound[nodeId][merPos]
		boolean[][] merFound = new boolean[session.ARTree.getNodeCount()][sk.getMerCount()];
		// loop on words
		byte[] qw = null;
		while ((qw = sk.getNextByteWord()) != null) {
			// get Pairs associated to this word
			Char2FloatMap.FastEntrySet allPairs = null;
			// TODO: what is this DNAStatesShifted
			if (session.states instanceof DNAStatesShifted) {
				allPairs = session.hash.getPairsOfTopPosition2(session.states.compressMer(qw));
			} else {
				allPairs = session.hash.getPairsOfTopPosition2(qw);
			}

			// word is not present in hash
			if (allPairs == null) {
				queryKmerCount++;
				continue;
			}
			queryKmerMatchingDB++;

			// stream version, 5-10% faster than allPairs.fastIterator()
			allPairs.stream().forEach((Char2FloatMap.Entry entry) -> {
				int nodeId = entry.getCharKey();
				// we will score only encountered nodes, originalNode registered
				// at 1st encouter
				if (nodeOccurences[nodeId] == 0) {
					selectedNodes.add(nodeId);
				}
				// count # times originalNode encountered
				nodeOccurences[nodeId] += 1;
				// score associated to originalNode x for current read
				nodeScores[nodeId] += entry.getFloatValue();
			});

			queryKmerCount++;
		}

		if (merStats) {
			try {
				PlacementProcess.merStats(session, merFound, fasta, bwMerStats);
			} catch (IOException e) {
				e.printStackTrace();
			}
		}
	}
}
