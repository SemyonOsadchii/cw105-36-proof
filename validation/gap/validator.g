#############################################################################
##
##  Clean-room GAP validation for the length-35 circulant sequence problem.
##  This file is self-contained and uses only GAP core finite-field,
##  polynomial, matrix, and integer arithmetic.
##
#############################################################################

TASK_ID := "cgw35-gap-publication-v1";
N := 35;
TARGET_WEIGHT := 12;

AssertTrue := function(condition, message)
    if not condition then
        Error(Concatenation("ASSERTION FAILED: ", message));
    fi;
end;

Hex35 := function(value)
    local digits, chars, x, i, r;
    digits := "0123456789abcdef";
    chars := [];
    x := value;
    for i in [1..9] do
        r := x mod 16;
        Add(chars, digits[r + 1]);
        x := QuoInt(x, 16);
    od;
    AssertTrue(x = 0, "35-bit mask exceeded nine hexadecimal digits");
    return Reversed(String(chars));
end;

JsonBool := function(value)
    if value then
        return "true";
    fi;
    return "false";
end;

JsonIntList := function(values)
    return Concatenation("[", JoinStringsWithSeparator(List(values, String), ","), "]");
end;

PolynomialFromCoefficients := function(field, x, coefficients)
    local result, i;
    result := Zero(x);
    for i in [1..Length(coefficients)] do
        result := result + coefficients[i] * x^(i - 1);
    od;
    return result;
end;

SharpPolynomial := function(field, x, polynomial, conjugationPower)
    local coefficients, reversed, leading;
    coefficients := CoefficientsOfUnivariatePolynomial(polynomial);
    reversed := List(Reversed(coefficients), c -> c^conjugationPower);
    leading := reversed[Length(reversed)];
    AssertTrue(leading <> Zero(field), "reciprocal polynomial has zero leading coefficient");
    reversed := List(reversed, c -> c / leading);
    return PolynomialFromCoefficients(field, x, reversed);
end;

ProductPolynomials := function(x, polynomials)
    local result, f;
    result := One(x);
    for f in polynomials do
        result := result * f;
    od;
    return result;
end;

FactorOrbits := function(field, x, factors, conjugationPower)
    local images, i, sharp, position, seen, orbits, orbit;
    images := [];
    for i in [1..Length(factors)] do
        sharp := SharpPolynomial(field, x, factors[i], conjugationPower);
        position := Position(factors, sharp);
        AssertTrue(position <> fail, "sharp image missing from factor list");
        Add(images, position);
    od;
    for i in [1..Length(images)] do
        AssertTrue(images[images[i]] = i, "factor sharp action is not an involution");
    od;
    seen := ListWithIdenticalEntries(Length(factors), false);
    orbits := [];
    for i in [1..Length(factors)] do
        if not seen[i] then
            if images[i] = i then
                orbit := [i];
                seen[i] := true;
            else
                orbit := Set([i, images[i]]);
                seen[orbit[1]] := true;
                seen[orbit[2]] := true;
            fi;
            Add(orbits, orbit);
        fi;
    od;
    return rec(images := images, orbits := orbits);
end;

MinimalBranchGenerators := function(x, factors, factorOrbits)
    local fixed, paired, orbit, branchCount, branches, choice, chosen, j, g;
    fixed := Filtered(factorOrbits, o -> Length(o) = 1);
    paired := Filtered(factorOrbits, o -> Length(o) = 2);
    branchCount := 2^Length(paired);
    branches := [];
    for choice in [0..branchCount - 1] do
        chosen := [];
        for orbit in fixed do
            Add(chosen, orbit[1]);
        od;
        for j in [1..Length(paired)] do
            Add(chosen, paired[j][1 + ((QuoInt(choice, 2^(j - 1))) mod 2)]);
        od;
        Sort(chosen);
        g := ProductPolynomials(x, List(chosen, i -> factors[i]));
        Add(branches, rec(
            choice := choice,
            factorIndices := chosen,
            generator := g,
            dimension := N - Degree(g)
        ));
    od;
    return rec(fixed := fixed, paired := paired, branches := branches);
end;

PolynomialCoefficientCodes3 := function(polynomial)
    return List(CoefficientsOfUnivariatePolynomial(polynomial), Int);
end;

Code4 := function(value, field, alpha)
    if value = Zero(field) then
        return 0;
    elif value = One(field) then
        return 1;
    elif value = alpha then
        return 2;
    elif value = alpha^2 then
        return 3;
    fi;
    Error("unexpected GF(4) element");
end;

PolynomialCoefficientCodes4 := function(polynomial, field, alpha)
    return List(CoefficientsOfUnivariatePolynomial(polynomial),
                c -> Code4(c, field, alpha));
end;

MaskOfWord := function(word, zero)
    local mask, q;
    mask := 0;
    for q in [0..N - 1] do
        if word[q + 1] <> zero then
            mask := mask + 2^q;
        fi;
    od;
    return mask;
end;

SupportOfMask := function(mask)
    return Filtered([0..N - 1], q -> (QuoInt(mask, 2^q) mod 2) = 1);
end;

ShiftedGeneratorRows := function(generator, field)
    local coefficients, degree, dimension, rows, j, row, i;
    coefficients := CoefficientsOfUnivariatePolynomial(generator);
    degree := Degree(generator);
    dimension := N - degree;
    rows := [];
    for j in [0..dimension - 1] do
        row := ListWithIdenticalEntries(N, Zero(field));
        for i in [0..degree] do
            row[j + i + 1] := coefficients[i + 1];
        od;
        Add(rows, row);
    od;
    return rows;
end;

WordFromGeneratorMap := function(generatorMap, message, field)
    local word, j, q, scalar;
    word := ListWithIdenticalEntries(N, Zero(field));
    for j in [1..Length(message)] do
        scalar := message[j];
        if scalar <> Zero(field) then
            for q in [1..N] do
                word[q] := word[q] + scalar * generatorMap[j][q];
            od;
        fi;
    od;
    return word;
end;

Char3CorrelationsZero := function(word, field)
    local t, q, sum;
    for t in [0..N - 1] do
        sum := Zero(field);
        for q in [0..N - 1] do
            sum := sum + word[q + 1] * word[((q + t) mod N) + 1];
        od;
        if sum <> Zero(field) then
            return false;
        fi;
    od;
    return true;
end;

SignCodeOfNormalizedWord := function(word, field)
    local support, firstValue, scale, code, j;
    support := Filtered([0..N - 1], q -> word[q + 1] <> Zero(field));
    AssertTrue(Length(support) = TARGET_WEIGHT, "sign encoding requires weight 12");
    firstValue := word[support[1] + 1];
    scale := One(field) / firstValue;
    code := 0;
    for j in [1..Length(support)] do
        if scale * word[support[j] + 1] = -One(field) then
            code := code + 2^(j - 1);
        else
            AssertTrue(scale * word[support[j] + 1] = One(field),
                       "normalized GF(3) entry is not a sign");
        fi;
    od;
    AssertTrue((code mod 2) = 0, "first normalized sign is not +1");
    return code;
end;

Char3WordFromMaskAndSignCode := function(mask, signCode, field)
    local word, support, j;
    word := ListWithIdenticalEntries(N, Zero(field));
    support := SupportOfMask(mask);
    for j in [1..Length(support)] do
        if (QuoInt(signCode, 2^(j - 1)) mod 2) = 0 then
            word[support[j] + 1] := One(field);
        else
            word[support[j] + 1] := -One(field);
        fi;
    od;
    return word;
end;

EnumerateChar3Branch := function(branch, field)
    local rows, dimension, digits, word, rawKeys, weightCount, total,
          step, carry, j, q, weight, mask, signCode, key;
    rows := ShiftedGeneratorRows(branch.generator, field);
    dimension := branch.dimension;
    AssertTrue(Length(rows) = dimension, "GF(3) generator row count mismatch");
    digits := ListWithIdenticalEntries(dimension, 0);
    word := ListWithIdenticalEntries(N, Zero(field));
    rawKeys := [];
    weightCount := 0;
    total := 3^dimension;
    for step in [1..total - 1] do
        carry := true;
        j := 1;
        while carry do
            for q in [1..N] do
                word[q] := word[q] + rows[j][q];
            od;
            if digits[j] = 2 then
                digits[j] := 0;
                j := j + 1;
            else
                digits[j] := digits[j] + 1;
                carry := false;
            fi;
        od;
        weight := Number(word, c -> c <> Zero(field));
        if weight = TARGET_WEIGHT then
            weightCount := weightCount + 1;
            mask := MaskOfWord(word, Zero(field));
            signCode := SignCodeOfNormalizedWord(word, field);
            key := mask * 2^TARGET_WEIGHT + signCode;
            Add(rawKeys, key);
        fi;
    od;
    return rec(rawKeys := rawKeys, weightCount := weightCount, enumerated := total);
end;

CompressNormalizedChar3Keys := function(rawKeys)
    local normalized, i, j;
    Sort(rawKeys);
    normalized := [];
    i := 1;
    while i <= Length(rawKeys) do
        j := i + 1;
        while j <= Length(rawKeys) and rawKeys[j] = rawKeys[i] do
            j := j + 1;
        od;
        AssertTrue(j - i = 2, "normalized GF(3) assignment did not occur as exactly {v,-v}");
        Add(normalized, rawKeys[i]);
        i := j;
    od;
    return normalized;
end;

WriteSupportMasks := function(filename, masks)
    local mask, stream;
    stream := OutputTextFile(filename, false);
    SetPrintFormattingStatus(stream, false);
    for mask in masks do
        PrintTo(stream, Hex35(mask), "\n");
    od;
    CloseStream(stream);
end;

AffineImageMask := function(mask, multiplier, translation)
    local image, q, target;
    image := 0;
    for q in [0..N - 1] do
        if (QuoInt(mask, 2^q) mod 2) = 1 then
            target := (multiplier * q + translation) mod N;
            image := image + 2^target;
        fi;
    od;
    return image;
end;

AffineSupportOrbits := function(masks, normalizedCounts)
    local units, marked, orbits, i, images, a, b, image, position,
          multiplicities, p;
    units := Filtered([0..N - 1], a -> GcdInt(a, N) = 1);
    AssertTrue(Length(units) = 24, "unit group modulo 35 was not derived with size 24");
    marked := ListWithIdenticalEntries(Length(masks), false);
    orbits := [];
    for i in [1..Length(masks)] do
        if not marked[i] then
            images := [];
            for a in units do
                for b in [0..N - 1] do
                    Add(images, AffineImageMask(masks[i], a, b));
                od;
            od;
            images := Set(images);
            multiplicities := [];
            for image in images do
                position := PositionSorted(masks, image);
                AssertTrue(position <= Length(masks) and masks[position] = image,
                           "affine image missing from frozen support family");
                marked[position] := true;
                Add(multiplicities, normalizedCounts[position]);
            od;
            multiplicities := Set(multiplicities);
            AssertTrue(Length(multiplicities) = 1,
                       "GF(3) assignment multiplicity is not affine invariant");
            Add(orbits, rec(
                representative := masks[i],
                size := Length(images),
                normalizedAssignmentsPerSupport := multiplicities[1],
                wordsPerSupport := 2 * multiplicities[1]
            ));
        fi;
    od;
    AssertTrue(ForAll(marked, x -> x), "affine orbit decomposition did not cover all supports");
    AssertTrue(Sum(List(orbits, o -> o.size)) = Length(masks),
               "affine orbit sizes do not sum to support count");
    return orbits;
end;

RightNullspaceBasis := function(matrix, ncols, field)
    local a, nrows, pivotColumns, row, col, pivot, i, factor,
          freeColumns, basis, free, vector, r;
    a := List(matrix, ShallowCopy);
    nrows := Length(a);
    pivotColumns := [];
    row := 1;
    for col in [1..ncols] do
        pivot := fail;
        for i in [row..nrows] do
            if a[i][col] <> Zero(field) then
                pivot := i;
                break;
            fi;
        od;
        if pivot <> fail then
            if pivot <> row then
                r := a[row];
                a[row] := a[pivot];
                a[pivot] := r;
            fi;
            factor := a[row][col];
            a[row] := List(a[row], x -> x / factor);
            for i in [1..nrows] do
                if i <> row and a[i][col] <> Zero(field) then
                    factor := a[i][col];
                    a[i] := a[i] - factor * a[row];
                fi;
            od;
            Add(pivotColumns, col);
            row := row + 1;
            if row > nrows then
                break;
            fi;
        fi;
    od;
    freeColumns := Difference([1..ncols], pivotColumns);
    basis := [];
    for free in freeColumns do
        vector := ListWithIdenticalEntries(ncols, Zero(field));
        vector[free] := One(field);
        for i in [1..Length(pivotColumns)] do
            vector[pivotColumns[i]] := -a[i][free];
        od;
        Add(basis, vector);
    od;
    return basis;
end;

GeneratorCoordinateMatrix := function(generator, field)
    local coefficients, degree, dimension, matrix, q, row, j, index;
    coefficients := CoefficientsOfUnivariatePolynomial(generator);
    degree := Degree(generator);
    dimension := N - degree;
    matrix := [];
    for q in [0..N - 1] do
        row := ListWithIdenticalEntries(dimension, Zero(field));
        for j in [0..dimension - 1] do
            index := q - j;
            if index >= 0 and index <= degree then
                row[j + 1] := coefficients[index + 1];
            fi;
        od;
        Add(matrix, row);
    od;
    return matrix;
end;

MessageToWordByCoordinateMatrix := function(coordinateMatrix, message, field)
    local word, q, j, sum;
    word := [];
    for q in [1..N] do
        sum := Zero(field);
        for j in [1..Length(message)] do
            sum := sum + coordinateMatrix[q][j] * message[j];
        od;
        Add(word, sum);
    od;
    return word;
end;

Char4CorrelationsZero := function(word, field)
    local t, q, sum;
    for t in [0..N - 1] do
        sum := Zero(field);
        for q in [0..N - 1] do
            sum := sum + word[q + 1] * word[((q + t) mod N) + 1]^2;
        od;
        if sum <> Zero(field) then
            return false;
        fi;
    od;
    return true;
end;

ExponentCodeOfNormalizedWord := function(word, support, field, alpha)
    local scale, code, j, value, digit;
    scale := One(field) / word[support[1] + 1];
    code := 0;
    for j in [1..Length(support)] do
        value := scale * word[support[j] + 1];
        if value = One(field) then
            digit := 0;
        elif value = alpha then
            digit := 1;
        elif value = alpha^2 then
            digit := 2;
        else
            Error("normalized GF(4) value is not in its multiplicative group");
        fi;
        code := code + digit * 3^(j - 1);
    od;
    AssertTrue((code mod 3) = 0, "first normalized GF(4) exponent is not zero");
    return code;
end;

Char4WordFromMaskAndExponentCode := function(mask, exponentCode, field, alpha)
    local support, word, value, j, digit;
    support := SupportOfMask(mask);
    word := ListWithIdenticalEntries(N, Zero(field));
    value := exponentCode;
    for j in [1..Length(support)] do
        digit := value mod 3;
        value := QuoInt(value, 3);
        word[support[j] + 1] := alpha^digit;
    od;
    return word;
end;

EnumerateChar4SupportBranch := function(mask, branch, coordinateMatrix, field, alpha)
    local support, outside, constraints, basis, dimension, fieldElements,
          rawKeys, vectorsExamined, number, temp, coefficients, i, j,
          message, word, wordMask, exponentCode, key, scalar;
    support := SupportOfMask(mask);
    outside := Difference([0..N - 1], support);
    constraints := List(outside, q -> coordinateMatrix[q + 1]);
    basis := RightNullspaceBasis(constraints, branch.dimension, field);
    dimension := Length(basis);
    fieldElements := [Zero(field), One(field), alpha, alpha^2];
    rawKeys := [];
    vectorsExamined := 0;
    for number in [1..4^dimension - 1] do
        vectorsExamined := vectorsExamined + 1;
        temp := number;
        coefficients := [];
        for i in [1..dimension] do
            Add(coefficients, fieldElements[(temp mod 4) + 1]);
            temp := QuoInt(temp, 4);
        od;
        message := ListWithIdenticalEntries(branch.dimension, Zero(field));
        for i in [1..dimension] do
            scalar := coefficients[i];
            if scalar <> Zero(field) then
                for j in [1..branch.dimension] do
                    message[j] := message[j] + scalar * basis[i][j];
                od;
            fi;
        od;
        word := MessageToWordByCoordinateMatrix(coordinateMatrix, message, field);
        wordMask := MaskOfWord(word, Zero(field));
        AssertTrue(IsSubsetSet(support, SupportOfMask(wordMask)),
                   "outside-support GF(4) constraint failed");
        if wordMask = mask then
            AssertTrue(Char4CorrelationsZero(word, field),
                       "branch-derived GF(4) word failed full Hermitian correlations");
            exponentCode := ExponentCodeOfNormalizedWord(word, support, field, alpha);
            key := mask * 3^TARGET_WEIGHT + exponentCode;
            Add(rawKeys, key);
        fi;
    od;
    return rec(
        kernelDimension := dimension,
        vectorsExamined := vectorsExamined,
        rawKeys := rawKeys
    );
end;

PairAdd := function(x, y)
    return [x[1] + y[1], x[2] + y[2]];
end;

PairMultiply := function(x, y)
    return [
        x[1] * y[1] - x[2] * y[2],
        x[1] * y[2] + x[2] * y[1] - x[2] * y[2]
    ];
end;

PairConjugate := function(x)
    return [x[1] - x[2], -x[2]];
end;

IsZeroPair := function(x)
    return x[1] = 0 and x[2] = 0;
end;

OriginalWordFromCodes := function(mask, signCode, exponentCode)
    local units, support, word, signTemp, exponentTemp, j, sign, exponent, index;
    units := [[1,0], [0,1], [-1,-1], [-1,0], [0,-1], [1,1]];
    support := SupportOfMask(mask);
    word := List([1..N], i -> [0,0]);
    signTemp := signCode;
    exponentTemp := exponentCode;
    for j in [1..Length(support)] do
        sign := signTemp mod 2;
        signTemp := QuoInt(signTemp, 2);
        exponent := exponentTemp mod 3;
        exponentTemp := QuoInt(exponentTemp, 3);
        index := exponent + 1 + 3 * sign;
        word[support[j] + 1] := ShallowCopy(units[index]);
    od;
    return word;
end;

OriginalDenseCandidateCheck := function(word)
    local t, q, sum, expected;
    for t in [0..N - 1] do
        sum := [0,0];
        for q in [0..N - 1] do
            sum := PairAdd(sum, PairMultiply(
                word[q + 1],
                PairConjugate(word[((q + t) mod N) + 1])
            ));
        od;
        if t = 0 then
            expected := [TARGET_WEIGHT,0];
        else
            expected := [0,0];
        fi;
        if sum <> expected then
            return false;
        fi;
    od;
    return true;
end;

VerifyDenseIndependent := function(word)
    local correlations, shift, left, right, accumulator;
    correlations := [];
    for shift in [0..34] do
        accumulator := [0,0];
        for left in [1..35] do
            right := ((left - 1 + shift) mod 35) + 1;
            accumulator := PairAdd(accumulator,
                PairMultiply(word[left], PairConjugate(word[right])));
        od;
        Add(correlations, accumulator);
    od;
    return correlations[1] = [12,0]
           and ForAll(correlations{[2..35]}, IsZeroPair);
end;

VerifySparseOrderedDifferences := function(word)
    local support, bins, i, j, shift, contribution;
    support := Filtered([0..34], q -> not IsZeroPair(word[q + 1]));
    bins := List([1..35], k -> [0,0]);
    for i in support do
        for j in support do
            shift := (j - i) mod 35;
            contribution := PairMultiply(word[i + 1], PairConjugate(word[j + 1]));
            bins[shift + 1] := PairAdd(bins[shift + 1], contribution);
        od;
    od;
    return bins[1] = [12,0] and ForAll(bins{[2..35]}, IsZeroPair);
end;

VerifyCRTGrid := function(word)
    local grid, q, r5, r7, s5, s7, a, b, accumulator, expected;
    grid := List([1..5], i -> List([1..7], j -> [0,0]));
    for q in [0..34] do
        r5 := q mod 5;
        r7 := q mod 7;
        grid[r5 + 1][r7 + 1] := ShallowCopy(word[q + 1]);
    od;
    for s5 in [0..4] do
        for s7 in [0..6] do
            accumulator := [0,0];
            for a in [0..4] do
                for b in [0..6] do
                    accumulator := PairAdd(accumulator, PairMultiply(
                        grid[a + 1][b + 1],
                        PairConjugate(grid[((a + s5) mod 5) + 1]
                                           [((b + s7) mod 7) + 1])
                    ));
                od;
            od;
            if s5 = 0 and s7 = 0 then
                expected := [12,0];
            else
                expected := [0,0];
            fi;
            if accumulator <> expected then
                return false;
            fi;
        od;
    od;
    return true;
end;

VerifyQuotientGroupRing := function(word)
    local sharp, q, product, i, j, degree, expected;
    sharp := List([1..35], k -> [0,0]);
    for q in [0..34] do
        sharp[((-q) mod 35) + 1] := PairConjugate(word[q + 1]);
    od;
    product := List([1..35], k -> [0,0]);
    for i in [0..34] do
        for j in [0..34] do
            degree := (i + j) mod 35;
            product[degree + 1] := PairAdd(product[degree + 1],
                PairMultiply(word[i + 1], sharp[j + 1]));
        od;
    od;
    for degree in [0..34] do
        if degree = 0 then
            expected := [12,0];
        else
            expected := [0,0];
        fi;
        if product[degree + 1] <> expected then
            return false;
        fi;
    od;
    return true;
end;

VerifyPhases := function(word)
    local units, nonzero, coefficient;
    units := [[1,0], [0,1], [-1,-1], [-1,0], [0,-1], [1,1]];
    nonzero := 0;
    for coefficient in word do
        if not IsZeroPair(coefficient) then
            if Position(units, coefficient) = fail then
                return false;
            fi;
            nonzero := nonzero + 1;
        fi;
    od;
    return nonzero = 12;
end;

CoefficientLabels := function(word)
    local units, labels, result, coefficient, position;
    units := [[1,0], [0,1], [-1,-1], [-1,0], [0,-1], [1,1]];
    labels := ["1", "omega", "omega^2", "-1", "-omega", "-omega^2"];
    result := [];
    for coefficient in word do
        if IsZeroPair(coefficient) then
            Add(result, "0");
        else
            position := Position(units, coefficient);
            AssertTrue(position <> fail, "cannot label non-unit coefficient");
            Add(result, labels[position]);
        fi;
    od;
    return result;
end;

FactorJson3 := function(factors)
    local objects, i;
    objects := [];
    for i in [1..Length(factors)] do
        Add(objects, Concatenation(
            "{\"index\":", String(i),
            ",\"degree\":", String(Degree(factors[i])),
            ",\"coefficients\":", JsonIntList(PolynomialCoefficientCodes3(factors[i])),
            "}"
        ));
    od;
    return Concatenation("[", JoinStringsWithSeparator(objects, ","), "]");
end;

FactorJson4 := function(factors, field, alpha)
    local objects, i;
    objects := [];
    for i in [1..Length(factors)] do
        Add(objects, Concatenation(
            "{\"index\":", String(i),
            ",\"degree\":", String(Degree(factors[i])),
            ",\"coefficients\":", JsonIntList(
                PolynomialCoefficientCodes4(factors[i], field, alpha)),
            "}"
        ));
    od;
    return Concatenation("[", JoinStringsWithSeparator(objects, ","), "]");
end;

OrbitJson := function(orbits)
    local objects, orbit;
    objects := [];
    for orbit in orbits do
        Add(objects, JsonIntList(orbit));
    od;
    return Concatenation("[", JoinStringsWithSeparator(objects, ","), "]");
end;

BranchJson := function(branches, enumerationCounts)
    local objects, i, extra;
    objects := [];
    for i in [1..Length(branches)] do
        if enumerationCounts = fail then
            extra := "";
        else
            extra := Concatenation(
                ",\"enumerated_words\":", String(enumerationCounts[i].enumerated),
                ",\"weight_12_words\":", String(enumerationCounts[i].weightCount)
            );
        fi;
        Add(objects, Concatenation(
            "{\"choice\":", String(branches[i].choice),
            ",\"factor_indices\":", JsonIntList(branches[i].factorIndices),
            ",\"generator_degree\":", String(Degree(branches[i].generator)),
            ",\"dimension\":", String(branches[i].dimension),
            extra, "}"
        ));
    od;
    return Concatenation("[", JoinStringsWithSeparator(objects, ","), "]");
end;

AffineOrbitJson := function(affineOrbits)
    local objects, orbit;
    objects := [];
    for orbit in affineOrbits do
        Add(objects, Concatenation(
            "{\"representative\":\"", Hex35(orbit.representative),
            "\",\"size\":", String(orbit.size),
            ",\"normalized_assignments_per_support\":",
            String(orbit.normalizedAssignmentsPerSupport),
            ",\"words_per_support\":", String(orbit.wordsPerSupport), "}"
        ));
    od;
    return Concatenation("[", JoinStringsWithSeparator(objects, ","), "]");
end;

WriteChar3Frozen := function(filename, factors, factorData, branchData,
                             branchCounts, supportMasks, normalizedKeys,
                             normalizedCounts, affineOrbits, directChecks)
    local stream;
    stream := OutputTextFile(filename, false);
    SetPrintFormattingStatus(stream, false);
    PrintTo(stream,
        "{\n",
        "\"task_id\":\"", TASK_ID, "\",\n",
        "\"frozen_stage\":\"characteristic_3_complete_support_family\",\n",
        "\"field\":\"GF(3)\",\n",
        "\"support_preservation_checked\":true,\n",
        "\"modulus_degree\":35,\n",
        "\"factorization\":", FactorJson3(factors), ",\n",
        "\"factor_sharp_images\":", JsonIntList(factorData.images), ",\n",
        "\"factor_sharp_orbits\":", OrbitJson(factorData.orbits), ",\n",
        "\"irreducibility_checked\":true,\n",
        "\"square_free_checked\":true,\n",
        "\"reconstruction_checked\":true,\n",
        "\"branch_lemma_checked\":true,\n",
        "\"branches\":", BranchJson(branchData.branches, branchCounts), ",\n",
        "\"total_enumerated_words\":",
        String(Sum(List(branchCounts, x -> x.enumerated))), ",\n",
        "\"total_weight_12_words\":",
        String(Sum(List(branchCounts, x -> x.weightCount))), ",\n",
        "\"normalized_weight_12_assignments\":", String(Length(normalizedKeys)), ",\n",
        "\"support_count\":", String(Length(supportMasks)), ",\n",
        "\"support_masks_file\":\"support-masks.txt\",\n",
        "\"support_masks_sorted_unique_fixed_width\":true,\n",
        "\"word_multiplicities_by_support\":",
        JsonIntList(List(normalizedCounts, x -> 2 * x)), ",\n",
        "\"direct_full_correlation_checks\":", String(directChecks), ",\n",
        "\"affine_group_order\":840,\n",
        "\"affine_orbits\":", AffineOrbitJson(affineOrbits), ",\n",
        "\"complete\":true\n",
        "}\n"
    );
    CloseStream(stream);
end;

WriteChar4Frozen := function(filename, factors, factorData, branchData,
                             supportCount, processedCases, dimensionHistogram,
                             vectorsExamined, rawOccurrenceCount, uniqueKeys,
                             survivorMasks, survivorCounts, directChecks)
    local survivorObjects, i, stream;
    survivorObjects := [];
    for i in [1..Length(survivorMasks)] do
        Add(survivorObjects, Concatenation(
            "{\"mask\":\"", Hex35(survivorMasks[i]),
            "\",\"normalized_assignments\":", String(survivorCounts[i]), "}"
        ));
    od;
    stream := OutputTextFile(filename, false);
    SetPrintFormattingStatus(stream, false);
    PrintTo(stream,
        "{\n",
        "\"task_id\":\"", TASK_ID, "\",\n",
        "\"frozen_stage\":\"characteristic_2_complete_assignment_family\",\n",
        "\"field\":\"GF(4)\",\n",
        "\"coefficient_encoding\":\"0,1,alpha,alpha^2\",\n",
        "\"conjugation\":\"z_to_z_squared\",\n",
        "\"conjugation_checked\":true,\n",
        "\"support_preservation_checked\":true,\n",
        "\"modulus_degree\":35,\n",
        "\"factorization\":", FactorJson4(factors, GF(4), Z(4)), ",\n",
        "\"factor_sharp_images\":", JsonIntList(factorData.images), ",\n",
        "\"factor_sharp_orbits\":", OrbitJson(factorData.orbits), ",\n",
        "\"irreducibility_checked\":true,\n",
        "\"square_free_checked\":true,\n",
        "\"reconstruction_checked\":true,\n",
        "\"branch_lemma_checked\":true,\n",
        "\"branches\":", BranchJson(branchData.branches, fail), ",\n",
        "\"input_support_count\":", String(supportCount), ",\n",
        "\"support_branch_cases_expected\":",
        String(supportCount * Length(branchData.branches)), ",\n",
        "\"support_branch_cases_processed\":", String(processedCases), ",\n",
        "\"kernel_dimension_histogram_d0_through_d17\":",
        JsonIntList(dimensionHistogram), ",\n",
        "\"nonzero_kernel_vectors_examined\":", String(vectorsExamined), ",\n",
        "\"raw_normalized_branch_occurrences\":", String(rawOccurrenceCount), ",\n",
        "\"unique_normalized_assignments\":", String(Length(uniqueKeys)), ",\n",
        "\"surviving_support_count\":", String(Length(survivorMasks)), ",\n",
        "\"eliminated_support_count\":", String(supportCount - Length(survivorMasks)), ",\n",
        "\"survivors\":[", JoinStringsWithSeparator(survivorObjects, ","), "],\n",
        "\"selected_correlation_shifts\":",
        JsonIntList([0..34]), ",\n",
        "\"full_hermitian_correlation_checks\":", String(directChecks), ",\n",
        "\"normalization\":\"global_GF4_nonzero_scalar_sets_first_support_value_to_1\",\n",
        "\"normalization_orbit_size\":3,\n",
        "\"complete\":true\n",
        "}\n"
    );
    CloseStream(stream);
end;

#############################################################################
## Stage 0: exact reductions and support preservation.
#############################################################################

startRuntime := Runtime();
F3 := GF(3);
one3 := One(F3);
F4 := GF(4);
alpha := Z(4);
one4 := One(F4);
eisensteinUnits := [[1,0], [0,1], [-1,-1], [-1,0], [0,-1], [1,1]];
reductions3 := List(eisensteinUnits, u -> ((u[1] + u[2]) mod 3) * one3);
reductions4 := List(eisensteinUnits,
                    u -> (u[1] mod 2) * one4 + (u[2] mod 2) * alpha);
AssertTrue(ForAll(reductions3, z -> z <> Zero(F3)),
           "Eisenstein units did not remain nonzero modulo (1-omega)");
AssertTrue(ForAll(reductions4, z -> z <> Zero(F4)),
           "Eisenstein units did not remain nonzero modulo 2");
AssertTrue(alpha^2 + alpha + one4 = Zero(F4),
           "GF(4) primitive element does not satisfy alpha^2+alpha+1");
AssertTrue(ForAll(Filtered(Elements(F4), z -> z <> Zero(F4)),
                  z -> z * z^2 = one4),
           "GF(4) conjugation/norm check failed");

#############################################################################
## Stage 1: characteristic 3 factorization and complete support derivation.
#############################################################################

R3 := PolynomialRing(F3);
x3 := IndeterminatesOfPolynomialRing(R3)[1];
M3 := x3^N - One(R3);
factors3 := Factors(R3, M3);
AssertTrue(ProductPolynomials(x3, factors3) = M3,
           "GF(3) factor reconstruction failed");
AssertTrue(ForAll(factors3, IsIrreducibleRingElement),
           "GF(3) factor irreducibility check failed");
AssertTrue(Gcd(M3, Derivative(M3)) = One(M3),
           "GF(3) modulus is not square-free");
factorData3 := FactorOrbits(F3, x3, factors3, 1);
branchData3 := MinimalBranchGenerators(x3, factors3, factorData3.orbits);
for branch in branchData3.branches do
    AssertTrue((branch.generator *
                SharpPolynomial(F3, x3, branch.generator, 1)) mod M3 = Zero(M3),
               "GF(3) minimal branch does not imply vv#=0");
od;

branchCounts3 := [];
rawKeys3 := [];
for branch in branchData3.branches do
    Print("GF3_BRANCH_START choice=", branch.choice,
          " dimension=", branch.dimension, "\n");
    branchResult3 := EnumerateChar3Branch(branch, F3);
    Add(branchCounts3, branchResult3);
    Append(rawKeys3, branchResult3.rawKeys);
    Print("GF3_BRANCH_DONE choice=", branch.choice, "\n");
od;

normalizedKeys3 := CompressNormalizedChar3Keys(rawKeys3);
supportMasks := Set(List(normalizedKeys3,
                         key -> QuoInt(key, 2^TARGET_WEIGHT)));
normalizedCounts3 := ListWithIdenticalEntries(Length(supportMasks), 0);
signCodesBySupport := List([1..Length(supportMasks)], i -> []);
directChecks3 := 0;
for key in normalizedKeys3 do
    mask := QuoInt(key, 2^TARGET_WEIGHT);
    signCode := key mod 2^TARGET_WEIGHT;
    supportPosition := PositionSorted(supportMasks, mask);
    AssertTrue(supportPosition <= Length(supportMasks)
               and supportMasks[supportPosition] = mask,
               "normalized GF(3) key has unknown support");
    normalizedCounts3[supportPosition] := normalizedCounts3[supportPosition] + 1;
    Add(signCodesBySupport[supportPosition], signCode);
    word3 := Char3WordFromMaskAndSignCode(mask, signCode, F3);
    AssertTrue(Char3CorrelationsZero(word3, F3),
               "normalized GF(3) word failed direct correlations");
    directChecks3 := directChecks3 + 1;
od;

# This file is the first frozen mathematical output.  No external value has
# been read or compared before it is completely written.
WriteSupportMasks("support-masks.txt", supportMasks);
affineOrbits := AffineSupportOrbits(supportMasks, normalizedCounts3);
WriteChar3Frozen("characteristic3-frozen.json", factors3, factorData3,
                 branchData3, branchCounts3, supportMasks, normalizedKeys3,
                 normalizedCounts3, affineOrbits, directChecks3);
Print("GF3_FROZEN\n");

#############################################################################
## Stage 2: characteristic 2/GF(4), constrained branch subspaces per support.
#############################################################################

R4 := PolynomialRing(F4);
x4 := IndeterminatesOfPolynomialRing(R4)[1];
M4 := x4^N - One(R4);
factors4 := Factors(R4, M4);
AssertTrue(ProductPolynomials(x4, factors4) = M4,
           "GF(4) factor reconstruction failed");
AssertTrue(ForAll(factors4, IsIrreducibleRingElement),
           "GF(4) factor irreducibility check failed");
AssertTrue(Gcd(M4, Derivative(M4)) = One(M4),
           "GF(4) modulus is not square-free");
factorData4 := FactorOrbits(F4, x4, factors4, 2);
branchData4 := MinimalBranchGenerators(x4, factors4, factorData4.orbits);
for branch in branchData4.branches do
    AssertTrue((branch.generator *
                SharpPolynomial(F4, x4, branch.generator, 2)) mod M4 = Zero(M4),
               "GF(4) minimal branch does not imply ww#=0");
    branch.coordinateMatrix := GeneratorCoordinateMatrix(branch.generator, F4);
od;

rawKeys4 := [];
processedCases4 := 0;
dimensionHistogram4 := ListWithIdenticalEntries(18, 0);
vectorsExamined4 := 0;
for maskIndex in [1..Length(supportMasks)] do
    if (maskIndex mod 100) = 0 then
        Print("GF4_SUPPORT_PROGRESS ", maskIndex, "/", Length(supportMasks), "\n");
    fi;
    mask := supportMasks[maskIndex];
    for branch in branchData4.branches do
        branchResult4 := EnumerateChar4SupportBranch(
            mask, branch, branch.coordinateMatrix, F4, alpha);
        processedCases4 := processedCases4 + 1;
        dimensionHistogram4[branchResult4.kernelDimension + 1] :=
            dimensionHistogram4[branchResult4.kernelDimension + 1] + 1;
        vectorsExamined4 := vectorsExamined4 + branchResult4.vectorsExamined;
        Append(rawKeys4, branchResult4.rawKeys);
    od;
od;
AssertTrue(processedCases4 =
           Length(supportMasks) * Length(branchData4.branches),
           "GF(4) support/branch cases were not completely processed");

rawOccurrenceCount4 := Length(rawKeys4);
uniqueKeys4 := Set(rawKeys4);
directChecks4 := 0;
survivorMasks4 := [];
survivorCounts4 := [];
lastMask := fail;
for key in uniqueKeys4 do
    mask := QuoInt(key, 3^TARGET_WEIGHT);
    exponentCode := key mod 3^TARGET_WEIGHT;
    if lastMask = fail or mask <> lastMask then
        Add(survivorMasks4, mask);
        Add(survivorCounts4, 0);
        lastMask := mask;
    fi;
    survivorCounts4[Length(survivorCounts4)] :=
        survivorCounts4[Length(survivorCounts4)] + 1;
    word4 := Char4WordFromMaskAndExponentCode(mask, exponentCode, F4, alpha);
    AssertTrue(MaskOfWord(word4, Zero(F4)) = mask,
               "decoded GF(4) assignment changed support");
    AssertTrue(Char4CorrelationsZero(word4, F4),
               "decoded GF(4) assignment failed full correlations");
    directChecks4 := directChecks4 + 1;
od;

# This is the second frozen mathematical output.  Exact lifting starts only
# after the complete GF(4) assignment family has been written.
WriteChar4Frozen("characteristic2-frozen.json", factors4, factorData4,
                 branchData4, Length(supportMasks), processedCases4,
                 dimensionHistogram4, vectorsExamined4, rawOccurrenceCount4,
                 uniqueKeys4, survivorMasks4, survivorCounts4, directChecks4);
Print("GF4_FROZEN\n");

#############################################################################
## Stage 3: exact lifts and final conclusion.
#############################################################################

liftPairsTested := 0;
witness := fail;
for key in uniqueKeys4 do
    if witness <> fail then
        break;
    fi;
    mask := QuoInt(key, 3^TARGET_WEIGHT);
    exponentCode := key mod 3^TARGET_WEIGHT;
    supportPosition := PositionSorted(supportMasks, mask);
    AssertTrue(supportPosition <= Length(supportMasks)
               and supportMasks[supportPosition] = mask,
               "GF(4) survivor support absent from GF(3) family");
    for signCode in signCodesBySupport[supportPosition] do
        liftPairsTested := liftPairsTested + 1;
        candidate := OriginalWordFromCodes(mask, signCode, exponentCode);
        if OriginalDenseCandidateCheck(candidate) then
            witness := candidate;
            break;
        fi;
    od;
od;

if witness = fail then
    finalStatus := "nonexistent";
    denseVerified := false;
    sparseVerified := false;
    crtVerified := false;
    groupRingVerified := false;
    phaseVerified := false;
    witnessLabels := [];
    liftComplete := true;
else
    denseVerified := VerifyDenseIndependent(witness);
    sparseVerified := VerifySparseOrderedDifferences(witness);
    crtVerified := VerifyCRTGrid(witness);
    groupRingVerified := VerifyQuotientGroupRing(witness);
    phaseVerified := VerifyPhases(witness);
    AssertTrue(denseVerified and sparseVerified and crtVerified
               and groupRingVerified and phaseVerified,
               "original witness failed an independent verifier");
    witnessLabels := CoefficientLabels(witness);
    finalStatus := "exists";
    liftComplete := false;
fi;

cpuMilliseconds := Runtime() - startRuntime;
if witness = fail then
    witnessJson := "null";
else
    witnessJson := Concatenation(
        "{\"coefficients\":[\"",
        JoinStringsWithSeparator(witnessLabels, "\",\""),
        "\"],\"dense_correlations\":true,\"sparse_ordered_differences\":true,",
        "\"crt_grid\":true,\"quotient_group_ring\":true,\"phase_verifier\":true}"
    );
fi;

resultStream := OutputTextFile("result.json", false);
SetPrintFormattingStatus(resultStream, false);
PrintTo(resultStream,
    "{\n",
    "\"task_id\":\"", TASK_ID, "\",\n",
    "\"status\":\"", finalStatus, "\",\n",
    "\"scope\":\"circulant_sequence_only\",\n",
    "\"length\":35,\n",
    "\"target_weight\":12,\n",
    "\"exact_arithmetic_only\":true,\n",
    "\"external_expected_values_compared\":false,\n",
    "\"characteristic_3\":{\n",
    "  \"factorization\":", FactorJson3(factors3), ",\n",
    "  \"factor_sharp_images\":", JsonIntList(factorData3.images), ",\n",
    "  \"factor_sharp_orbits\":", OrbitJson(factorData3.orbits), ",\n",
    "  \"branches\":", BranchJson(branchData3.branches, branchCounts3), ",\n",
    "  \"support_count\":", String(Length(supportMasks)), ",\n",
    "  \"weight_12_word_count\":",
    String(Sum(List(branchCounts3, x -> x.weightCount))), ",\n",
    "  \"normalized_assignment_count\":", String(Length(normalizedKeys3)), ",\n",
    "  \"affine_orbit_count\":", String(Length(affineOrbits)), ",\n",
    "  \"affine_group_order\":840,\n",
    "  \"affine_orbits\":", AffineOrbitJson(affineOrbits), ",\n",
    "  \"word_multiplicities_by_support\":",
    JsonIntList(List(normalizedCounts3, x -> 2 * x)), ",\n",
    "  \"direct_full_correlation_checks\":", String(directChecks3), ",\n",
    "  \"complete\":true\n",
    "},\n",
    "\"characteristic_2\":{\n",
    "  \"factorization\":", FactorJson4(factors4, F4, alpha), ",\n",
    "  \"factor_sharp_images\":", JsonIntList(factorData4.images), ",\n",
    "  \"factor_sharp_orbits\":", OrbitJson(factorData4.orbits), ",\n",
    "  \"branches\":", BranchJson(branchData4.branches, fail), ",\n",
    "  \"support_branch_cases_expected\":",
    String(Length(supportMasks) * Length(branchData4.branches)), ",\n",
    "  \"support_branch_cases_processed\":", String(processedCases4), ",\n",
    "  \"kernel_dimension_histogram_d0_through_d17\":",
    JsonIntList(dimensionHistogram4), ",\n",
    "  \"nonzero_kernel_vectors_examined\":", String(vectorsExamined4), ",\n",
    "  \"raw_normalized_branch_occurrences\":", String(rawOccurrenceCount4), ",\n",
    "  \"unique_normalized_assignment_count\":", String(Length(uniqueKeys4)), ",\n",
    "  \"surviving_support_count\":", String(Length(survivorMasks4)), ",\n",
    "  \"eliminated_support_count\":",
    String(Length(supportMasks) - Length(survivorMasks4)), ",\n",
    "  \"full_correlation_checks\":", String(directChecks4), ",\n",
    "  \"support_branch_case_completeness_checked\":true,\n",
    "  \"complete\":true\n",
    "},\n",
    "\"exact_lifting\":{\n",
    "  \"normalized_sign_phase_pairs_tested\":", String(liftPairsTested), ",\n",
    "  \"exhaustion_complete\":", JsonBool(liftComplete), "\n",
    "},\n",
    "\"support_masks_file\":\"support-masks.txt\",\n",
    "\"characteristic_3_frozen_file\":\"characteristic3-frozen.json\",\n",
    "\"characteristic_2_frozen_file\":\"characteristic2-frozen.json\",\n",
    "\"witness\":", witnessJson, ",\n",
    "\"gap_cpu_milliseconds\":", String(cpuMilliseconds), ",\n",
    "\"completeness_certificate\":{\n",
    "  \"support_preservation_both_reductions\":true,\n",
    "  \"gap_factorizations_reconstructed_irreducible_square_free\":true,\n",
    "  \"factor_sharp_orbits_derived\":true,\n",
    "  \"all_minimal_cyclic_code_branches_covered\":true,\n",
    "  \"all_characteristic_3_branch_words_enumerated\":true,\n",
    "  \"all_frozen_supports_processed_in_characteristic_2\":true,\n",
    "  \"all_characteristic_2_kernel_subspaces_enumerated\":true,\n",
    "  \"all_surviving_assignments_checked_against_full_correlations\":true,\n",
    "  \"final_decision_is_complete\":true\n",
    "}\n",
    "}\n"
);
CloseStream(resultStream);

Print("FINAL_STATUS=", finalStatus, "\n");
Print("RESULT_JSON_WRITTEN\n");
QUIT_GAP(0);
