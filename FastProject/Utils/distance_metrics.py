import numpy as np;


def weighted_distance(pts, weights):
    # Un-normalized
    result = np.zeros((pts.shape[0], pts.shape[0]));
    waA = pts * weights;  # hadamard product
    result += np.dot(waA * pts, weights.T);
    result += result.T.copy();  # Copy is essential, otherwise memory is shared and things get weird
    result -= 2 * np.dot(waA, waA.T);
    np.fill_diagonal(result, 0);
    return result;


def weighted_distance_between(ptsA, weightsA, ptsB, weightsB):
    """
    Distances between the rows of A and the rows of B, weighted
    """
    result = np.zeros((ptsA.shape[0], ptsB.shape[0]));
    waA = ptsA * weightsA;
    waB = ptsB * weightsB;
    result += np.dot(waA * ptsA, weightsB.T);
    result += np.dot(weightsA, (waB * ptsB).T);
    result -= 2 * np.dot(waA, waB.T);
    return result;


def weighted_distance_between_slow(ptsA, weightsA, ptsB, weightsB):  # For debugging
    result = np.zeros((ptsA.shape[0], ptsB.shape[0]));
    for i in range(ptsA.shape[0]):
        for j in range(ptsB.shape[0]):
            result[i, j] = np.sum((ptsA[i, :] - ptsB[j, :])**2 * weightsA[i, :] * weightsB[j, :]);

    return result;


def weighted_distance_normalized(pts, weights):
    un_normalized = weighted_distance(pts, weights);
    un_normalized /= np.dot(weights, weights.T);
    normalized = un_normalized;  # Just a name change :)
    normalized[~np.isfinite(normalized)] = 0.0;  # For cases where all weights are 0
    return normalized;


def modeled_distance(pts, weights):
    wmean = np.sum(pts * weights, axis=0, keepdims=True) / np.sum(weights, axis=0, keepdims=True);
    wvar = np.sum((pts - wmean)**2 * weights, axis=0, keepdims=True) / np.sum(weights, axis=0, keepdims=True);

    mumat = wmean * np.ones(pts.shape);
    stdmat = np.sqrt(wvar) * np.ones(pts.shape);

    result = np.zeros((pts.shape[0], pts.shape[0]));

    # Both bad measurements
    result += 2 * np.dot(stdmat * (1 - weights), (stdmat * (1 - weights)).T);

    # One bad measurement.
    temp = weighted_distance_between(pts, weights, mumat, (1 - weights));  # One bad measurement
    temp += np.dot(stdmat * weights, (stdmat * (1 - weights)).T);  # Correction for one bad measurement case

    result += temp;
    result += temp.T;  # Multiply again to capture the symmetric scenario

    # Both good measurements
    result += weighted_distance(pts, weights);  # Case where both are good measurements

    return result;
