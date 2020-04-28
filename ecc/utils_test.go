package ecc

import (
	"fmt"
	"math/big"
	"testing"
)

func TestNafDecomposition(t *testing.T) {
	for i := int64(1); i < 17; i++ {
		exp := big.NewInt(i)
		var result [400]int8
		lExp := NafDecomposition(exp, result[:])
		dec := result[:lExp]

		// res := [5]int8{1, 0, -1, 0, 1}
		// for i, v := range dec {
		// 	if v != res[i] {
		// 		t.Error("Error in NafDecomposition")
		// 	}
		// }
		fmt.Println(exp, dec)

	}
}
