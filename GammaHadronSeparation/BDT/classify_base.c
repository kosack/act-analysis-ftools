
int main(int argc, char *argv[])
{

  char identifier[100];
  double msw,msl,xmax,xmaxerr;
  FILE *fp;
  void *attr[4];
  double rvals[2];
  double *ret,val;
  char line[256];
  char tok[256];

  ret = rvals;

  attr[0] = (void*) &msw;
  attr[1] = (void*) &msl;
  attr[2] = (void*) &xmax;
  attr[3] = (void*) &xmaxerr;
  
  if (!(fp = fopen( argv[1],"r" ))) {
    perror("open");
    printf("Couldn't open file %s\n", argv[1]);
    return(1);
  }
  
  while (fp && !feof(fp)) {

    if (fscanf(fp,"%s",line)<1)
      break;

    strcpy(identifier, strtok(line, ",;"));
    msw = atof(strtok(NULL, ",;"));
    msl = atof(strtok(NULL, ",;"));
    xmax = atof(strtok(NULL, ",;"));
    xmaxerr = atof(strtok(NULL, ",;"));

    val = predict( attr, ret );
    //    printf("%s (%lf %lf %lf %lf) %lf\n",identifier,msw,msl,xmax,xmaxerr,val);
    printf("%lf\n",val);
  }
  
  fclose(fp);

  return 0;
}
