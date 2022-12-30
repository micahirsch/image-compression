# Image Compression

The idea behind this project is to implement a barely simplified version of the image codification algorithm for `.jpg` archives, which is based on an essential property of the Fourier Transform: transforming signals from real life usually leads to low frequencies with little contribution from high frequencies. This notion introduces the following idea: given a signal we can consider it's Fourier Transform and discard higher frequencies. This allows us to store only a small number of frequencies, the compressed signal. To retrieve the signal you can complete the transformed vector with zeros and do the inverse transform.

This algorithm uses the Discrete Cosine Transform, mainly due to the fact that it converts Real numbers into Real numbers.

## Preparation

This simplification of the algorithm asumes the image dimensions are divisibles by 16. Therefore, the the algorithm has to firstly fill the image's margins with black to ensure that both dimensions are divisible by 16.

## Using the right format

The algorithm uses the YCbCr decomposition. Therefore, first we must convert the image to the correct format. We also separate the three channels to operate with them properly. Moreover, due to the human eye being muvh more sensitive to the luminosity than to the intensity of a color, we can reduce the matrices Cb and Cr by creating smaller matrices in which every pixel is the average of 4 neighbouring pixels in the original matrix. Finally, we'll center all of the matrices coefficients on 0. The values will range between -128 and 128.

The inverse process is also implemented, since it will be used later to retrieve the compressed image.

## Using the Discrete Cosine Transform

The following step consists in thinking of each matrix as a group of `8×8` submatrices and applying the discrete fourier transform to each.

Similarly to the previous step, the inverse process is also implemented, since it will be used later to retrieve the compressed image.

## Quantization

This is the most important compression step, in which we discard the highest frequencies in each `8×8` submatrix using a quantization matrix. This matrix can take different numbers and it's values change how the image is compressed. The quantization process consists in taking each submatrix and dividing it number to number by the quantization matrix rounding the result.

The result will have many 0s specially on the last rows and columns, which will help the compression later on.

## Compression 

Finally, we discard the zeros from the transformed submatrices. That is done by the following procedure:

1. Each submatrix is read in zig-zag and turned into a vector.

2. The resulting vector is compressed using the _Run Length Encoding_ method, which consists of indicating in a tuple the number of consecutive appearances of a number and the number itself. 

The vectors produced are all joined into a single vector.

## Joining the previous steps

The final algorithm is excecuted using the previously explained procedures and tested using different images and two different quantization matrices. The first quantization matrix is the most used matrix, meanwhile the second one was created by myself. 

The compression of an image using both matrices:
<p float="left", align= "middle">
  <img src="https://user-images.githubusercontent.com/83768210/210096327-83c7dc47-c167-4800-8c96-769a7dcbe8a2.png" width="500" />
  <img src="https://user-images.githubusercontent.com/83768210/210096377-08e04fd7-2a06-4fd6-91af-1910ea467920.png" width="500" /> 
</p>

Depending on which quantization matrix is used the image can preserve more or less quality.
