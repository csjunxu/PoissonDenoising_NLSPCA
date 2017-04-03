function plotimage(img)

    imagesc(img);
    axis image;
    caxis([0 255]);
   %caxis([0 max(max(img))]);
    colormap(gray(1024));
    axis off;
