            case 1
                % Instead of using a 2D gaussian kernel,
                % uses the fact that a gaussian filter can be separated in 1D gaussian kernels.
                A1 = imfilter(I1, reshape(Gx,[length(Gx) 1]),'symmetric');
                A1 = imfilter(shiftdim(A1,1), reshape(Gy,[length(Gy) 1]),'symmetric').';
                B1 = imfilter(I2, reshape(Gix,[length(Gix) 1]),'symmetric');
                B1 = imfilter(shiftdim(B1,1), reshape(Giy,[length(Giy) 1]),'symmetric').';
            case 2
                A1 = imfilter(I1, reshape(Gx,[length(Gx) 1]),'symmetric');
                A1 = imfilter(shiftdim(A1,1), reshape(Gdy,[length(Gdy) 1]),'symmetric').';
                B1 = imfilter(I2, reshape(Gix,[length(Gix) 1]),'symmetric');
                B1 = imfilter(shiftdim(B1,1), reshape(Gdiy,[length(Gdiy) 1]),'symmetric').';
            case 3
                A1 = imfilter(I1, reshape(Gdx,[length(Gdx) 1 1]),'symmetric');
                A1 = imfilter(shiftdim(A1,1), reshape(Gy,[length(Gy) 1 1]),'symmetric').';
                B1 = imfilter(I2, reshape(Gdix,[length(Gdix) 1]),'symmetric');
                B1 = imfilter(shiftdim(B1,1), reshape(Giy,[length(Giy) 1]),'symmetric').';
        end