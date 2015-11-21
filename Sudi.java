import java.io.*;
import java.util.*;

public class Sudi {
	private static final String start = "*";
	
	static Map<String,Map<String,Integer>> transitions; // prevPOS -> (nextPOS -> count)
	static TreeMap<String,Map<String,Integer>> emissions; // POS -> (word -> count)
	
	static Map<String,Map<String,Double>> transitionScores; // prevPOS -> (nextPOS -> score)
	static Map<String,Map<String,Double>> emissionScores; // currPOS -> (word -> score)
	
	/**
	 * Gets called by trainProbs and updates instance variables of transitions and emissions
	 * @param sentenceFile
	 * @param tagFile
	 * @throws Exception
	 */
	public static void trainCounts(String sentenceFile, String tagFile) throws Exception {		
		transitions = new TreeMap<String,Map<String,Integer>>();
		emissions = new TreeMap<String,Map<String,Integer>>();
		
		BufferedReader trainingSentences = new BufferedReader(new FileReader(sentenceFile));
		BufferedReader trainingTags = new BufferedReader(new FileReader(tagFile));
		
		String trainSent = trainingSentences.readLine();
		String trainTag = trainingTags.readLine();
		while (trainSent != null && trainTag != null) {  // Read each word and corresponding part of speech
			String[] trainWords = trainSent.toLowerCase().split(" ");
			String[] trainPOSs = trainTag.split(" ");
			
			for (int i = 0; i < trainWords.length && i < trainPOSs.length; i++) { // Scan up to the final period
				String prevPOS;
				if (i == 0) prevPOS = start;
				else prevPOS = trainPOSs[i-1];
				String nextPOS = trainPOSs[i];
				
				String currWord = trainWords[i];
				
				if (transitions.containsKey(prevPOS)) { // Check whether transitions already has prevPOS
					// Check whether the transition has not been made from the currPOS to the nextPOS
					if (!transitions.get(prevPOS).containsKey(nextPOS)) transitions.get(prevPOS).put(nextPOS, 1);
					// Otherwise update count of transition from currState to nextState
					else transitions.get(prevPOS).put(nextPOS, transitions.get(prevPOS).get(nextPOS) + 1);
				} else { // Otherwise add the currState and start count to the selected nextState
					Map<String,Integer> nextPOSToCount = new TreeMap<String,Integer>();
					nextPOSToCount.put(nextPOS, 1);
					transitions.put(prevPOS, nextPOSToCount);
				}
				
				if (emissions.containsKey(nextPOS)) { // Check whether emissions has the part of speech
					// Check whether the word has not been used in the given part of speech
					if (!emissions.get(nextPOS).containsKey(currWord)) emissions.get(nextPOS).put(currWord, 1);
					// Otherwise update count of the part of speech
					else emissions.get(nextPOS).put(currWord, emissions.get(nextPOS).get(currWord) + 1);
				} else { // Otherwise add the word and its part of speech
					Map<String,Integer> wordToCount = new TreeMap<String,Integer>();
					wordToCount.put(currWord, 1);
					emissions.put(nextPOS, wordToCount);
				}
			}
			// Go to next word and next part of speech
			trainSent = trainingSentences.readLine();
			trainTag = trainingTags.readLine();
		}
		trainingSentences.close();
		trainingTags.close();
	}
	
	/**
	 * Uses training files and trainCounts to teach Sudi how to recognize parts of speech within sentences
	 * @param sentenceFile
	 * @param tagFile
	 * @throws Exception
	 */
	public static void teachSudi(String sentenceFile, String tagFile) throws Exception {
		trainCounts(sentenceFile, tagFile); // Update the instance variables of transitions and emissions to calculate probabilities
		
		transitionScores = new TreeMap<String,Map<String,Double>>();
		emissionScores = new TreeMap<String,Map<String,Double>>();
		
		for (String prevPOS : transitions.keySet()) { // Scan through all states
			double total = 0;  // Calculate the total number of transitions from that state
			for (String nextPOS : transitions.get(prevPOS).keySet()) total += transitions.get(prevPOS).get(nextPOS); 
			
			// Calculate the probability of a transition from prevState to nextState
			Map<String,Double> prevToNextScore = new TreeMap<String,Double>();
			for (String nextPOS : transitions.get(prevPOS).keySet()) prevToNextScore.put(nextPOS, Math.log(transitions.get(prevPOS).get(nextPOS) / total));
			transitionScores.put(prevPOS, prevToNextScore);
		}

		for (String currPOS : emissions.keySet()) { // Scan through all observed emissions
			double total = 0;  // Calculate the number of ways that word is used as a part of speech
			for (String word : emissions.get(currPOS).keySet()) total += emissions.get(currPOS).get(word); 
			
			// Calculate the probability of the word being used as the given part of speech
			Map<String,Double> posToWordScore = new TreeMap<String,Double>();
			for (String word : emissions.get(currPOS).keySet()) posToWordScore.put(word, Math.log(emissions.get(currPOS).get(word) / total));
			emissionScores.put(currPOS, posToWordScore);
		}
	}
	
	/**
	 * Accepts a sentence and returns the parts of speech of the words in that sentence.
	 * Uses Viterbi's algorithm and depends on previous teachSudi
	 * @param sentence
	 * @return
	 */
	public static String viterbi(String sentence) {
		String[] input = sentence.toLowerCase().split(" ");	
		List<Map<String,String>> nextToPrevPOSs = new ArrayList<Map<String,String>>();  // Enables back-chaining at end
		
		Set<String> prevPOSs = new TreeSet<String>();  // Possible POSs of prevWord
		prevPOSs.add(start);
		
		Map<String,Double> prevScores = new TreeMap<String,Double>(); // Scores from prevWord
		prevScores.put(start, 0.0);
		
		for (int i = 0; i < input.length; i++) { // Loop through each word of input sentence
			Set<String> nextPOSs = new TreeSet<String>();  // Possible POSs of currWord
			Map<String,Double> nextScores = new TreeMap<String,Double>();  // Scores of currWord
			Map<String,String> nextToPrevPOS = new TreeMap<String,String>();    // Enters back-chaining at the end
			
			for (String prevPOS : prevPOSs) { // Loop through the prevPOSs of currWord
				if (transitionScores.containsKey(prevPOS)) {
					for (String nextPOS : transitionScores.get(prevPOS).keySet()) {
						nextPOSs.add(nextPOS);  // Add possiblePOS to nextPOSs
						
						// Calculate nextScore from currWord
						double emissionScore;
						// Check whether nextPOS emits currWord; if not, give it a high value -- unseen word penalty
						if (!emissionScores.get(nextPOS).containsKey(input[i])) emissionScore = -100;
						else emissionScore = emissionScores.get(nextPOS).get(input[i]);
						// Only considering possible nextPOSs
						double transitionScore = transitionScores.get(prevPOS).get(nextPOS);
						
						double nextScore = prevScores.get(prevPOS) + transitionScore + emissionScore;
						if (!nextScores.containsKey(nextPOS) || nextScore > nextScores.get(nextPOS)) {
							nextScores.put(nextPOS, nextScore);  // Updates possible nextScores to highest values
							nextToPrevPOS.put(nextPOS, prevPOS);  // Updates back-chain
						}
					}
				}
			}
			// Get ready for next currWord
			nextToPrevPOSs.add(nextToPrevPOS);
			prevPOSs = nextPOSs;
			prevScores = nextScores;
		}
		
		String partsOfSpeech = "";
		String currPOS = "";
		double bestScore = -500000.0; // High value that will be replaced
		for (String pos : prevScores.keySet()) {
			if (prevScores.get(pos) > bestScore) {
				bestScore = prevScores.get(pos);
				currPOS = pos; // endPOS has bestScore at end
			}
		}
		// Prepare POS statement to return
		for (int w = input.length - 1; w >= 0; w--) {
			partsOfSpeech = currPOS + " " + partsOfSpeech;
			currPOS = nextToPrevPOSs.get(w).get(currPOS);
		}
		return partsOfSpeech;
	}
	
	/**
	 * Prompts user to give Sudi a sentence to convert into 
	 * parts of speech and prints the conversion to the console.
	 */
	public static void giveSudi() {
		System.out.println("Hello, I'm Sudi. I convert the words of sentences into their parts of speech.");
		System.out.println("If you would like to quit, type 'quit'. Otherwise, give me a sentence:");
		Scanner userInput = new Scanner(System.in);
		String givenSentence = userInput.nextLine();
		while (!givenSentence.equals("quit")) {
			System.out.println(viterbi(givenSentence)+"\n--------------------");
			System.out.println("Give me another request:");
			givenSentence = userInput.nextLine();
		}
		System.out.println("Goodbye!");
	}
	
	/**
	 * Analyzes the Sudi's performance with two test files for the sentences and the tags
	 * @param sentenceFile
	 * @param tagFile
	 * @throws Exception
	 */
	public static void testSudi(String sentenceFile, String tagFile) throws Exception {
		BufferedReader testingSentences = new BufferedReader(new FileReader(sentenceFile));
		BufferedReader testingTags = new BufferedReader(new FileReader(tagFile));
		
		int right = 0;
		int wrong = 0;
		
		String testSent = testingSentences.readLine();
		String testTag = testingTags.readLine();
		while (testSent != null && testTag != null) {  // Read each word and corresponding part of speech
			String response = viterbi(testSent);
			
			String[] educatedGuess = response.split(" ");
			String[] answer = testTag.split(" ");
			
			for (int pos = 0; pos < educatedGuess.length && pos < answer.length; pos++) {
				if (educatedGuess[pos].equals(answer[pos])) right += 1;
				else wrong += 1;
			}
	
			testSent = testingSentences.readLine();
			testTag = testingTags.readLine();
		}
		double score = (double) right / (right + wrong);
		System.out.println("Test Score: " + score*100 + " %" + "\n---------");
		System.out.println("right: " + right);
		System.out.println("wrong: " + wrong);
		
		testingSentences.close();
		testingTags.close();
	}
	
	public static void main(String[] args) throws Exception {
		teachSudi("texts/simple-train-sentences.txt","texts/simple-train-tags.txt");
		testSudi("texts/simple-test-sentences.txt", "texts/simple-test-tags.txt");
		System.out.println("-------------------------------------");
		//System.out.println(viterbi("the dog saw trains in the night ."));
		//giveSudi();
	}
}