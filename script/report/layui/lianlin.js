function icon_append(){
    layui.use(['jquery','table'],function(){
        var $=layui.jquery;
        $('h1').prepend('<i class="layui-icon layui-icon-snowflake" style="font-size:30px;margin-right:15px;color:skyblue"></i>');
    })
    
}



function table_function(){
    var table=layui.table;
    table.init('demo', {height: 500,cellMinWidth:80, limit: 10,toolbar:true,page:true,defaultToolbar:['filter','exports'],});
    table.render();
    
}
    
function tab1(){
    var form=layui.form,$=layui.jquery;
    form.render('select');
        $('.tab1-quater').hide();
        var sample=$('#tab1').val();
        $('#tab1-'+sample).show();
        form.on('select(tab1)',function(data){
        $('.tab1-quater').hide();
            var sample=$('#tab1').val();
            $('#tab1-'+sample).show();
        })
    }

function tab2(){
    layui.use(['form','jquery'],function(){
    var form=layui.form,$=layui.jquery;
    form.render('select');
    $('.tab2-quater').hide();
    var sample=$('#tab2').val();
    $('#tab2-'+sample).show();
    form.on('select(tab2)',function(data){
    $('.tab2-quater').hide();
        var sample=$('#tab2').val();
        $('#tab2-'+sample).show();
    })
    
})
}

function tab3(){
    layui.use(['form','jquery'],function(){
    var form=layui.form,$=layui.jquery;
    form.render('select');
    $('.tab3-quater').hide();
    var sample=$('#tab3').val();
    $('#tab3-'+sample).show();
    form.on('select(tab3)',function(data){
    $('.tab3-quater').hide();
        var sample=$('#tab3').val();
        $('#tab3-'+sample).show();
    })
    
})
}


function shortvariant_stat(){
    var form=layui.form,$=layui.jquery;
    form.render();
    $('#snp').show();
    $('#indel').hide();
    form.on('radio(variant_stat)',function(data){
        if (data.value == 'snp'){
            $('#snp').show();
            $('#indel').hide();
        }
        else {
            $('#snp').hide();
            $('#indel').show();
        }
    })
}

function shortvariant_quality(){
    var form=layui.form,$=layui.jquery;
    form.render();
    $('.snp-quater').hide();
    $('.indel-quater').hide();
    var sample=$('#sample').val();
    $('#snp-'+sample).show();
    form.on('radio()',function(data){
        var variant=data.value;
        var sample=$('#sample').val();
        if (variant == 'snp'){
            $('#snp-'+sample).show();
            $('.indel-quater').hide();
        }
        else {
            $('#indel-'+sample).show();
            $('.snp-quater').hide();
        }
        $('#'+variant+'-'+sample).show();
    })

    form.on('select(sample)',function(data){
        var sample=data.value;
        //alert('重新选择SNP或者INDEL');
        var variant=$('#variant_select input:radio:checked').val();
        $('.snp-quater').hide();
        $('.indel-quater').hide();
        if (variant == 'snp'){
            $('#snp-'+sample).show();
        }
        else {
            $('#indel-'+sample).show();
        }
        $('#'+variant+'-'+sample).show();
        //alert(variant);
        form.on('radio(variant_check)',function(data){
            var variant=data.value;
            $('.snp-quater').hide();
            $('.indel-quater').hide();
            if (variant == 'snp'){
                $('#snp-'+sample).show();
            }
            else {
                $('#indel-'+sample).show();
            }

        })
    })
}

function shortvariant_func(){
    var form=layui.form,$=layui.jquery;
    form.render();
    $('#snp_func').show();
    $('#indel_func').hide();
    form.on('radio(variant_check)',function(data){
        if (data.value == 'snp'){
            $('#snp_func').show();
            $('#indel_func').hide();
        }
        else {
            $('#snp_func').hide();
            $('#indel_func').show();
        }
    })
}
