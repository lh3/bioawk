BEGIN{
    OFS=","	
}
{
    print $name,$comment,$seq
}
